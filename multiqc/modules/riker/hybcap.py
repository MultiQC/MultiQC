"""Parse riker `hybcap` outputs (hybcap-metrics.txt)."""

import logging
from typing import Dict, List

from multiqc import config
from multiqc.plots import linegraph, table
from multiqc.plots.table import TableConfig

from .util import read_tsv, to_float

log = logging.getLogger(__name__)

# Coverage thresholds emitted by riker as `frac_target_bases_{N}x`.
TARGET_COVERAGE_LEVELS = [1, 10, 20, 30, 50, 100, 250, 500, 1000]


def parse_reports(module):
    data_by_sample: Dict[str, Dict[str, float]] = {}
    panel_by_sample: Dict[str, str] = {}

    for f in module.find_log_files("riker/hybcap_metrics", filehandles=True):
        for row in read_tsv(f["f"]):
            sample = row.pop("sample", None)
            if not sample:
                continue
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="hybcap_metrics")

            panel = row.pop("panel_name", "")
            panel_by_sample[s_name] = panel

            data_by_sample[s_name] = {col: to_float(val) for col, val in row.items()}

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    riker_config = getattr(config, "riker_config", {})
    user_covs = riker_config.get("general_stats_target_coverage", [])
    if isinstance(user_covs, list) and user_covs:
        gs_covs: List[int] = [int(c) for c in user_covs if int(c) in TARGET_COVERAGE_LEVELS]
    else:
        gs_covs = [30]

    # Coloring rationale:
    #   - mean_target_coverage: experiment-dependent (panel + sequencing depth),
    #     so single-hue GnBu rather than red->green.
    #   - frac_selected_bases: higher is better; gradient over 50-100% so the
    #     practical range is differentiable.
    #   - frac_uncovered_targets: lower is better. Anything beyond ~10% is
    #     concerning, so cap the gradient at 10%.
    #   - hs_penalty_*: 1.0 is perfect efficiency. Cap gradient at 12 — typical
    #     well-behaved panels are 1-3, but values up to ~12 are reasonable in
    #     practice for deeper-coverage thresholds, so don't saturate too early.
    #   - at_dropout / gc_dropout: capped at 10 (Picard convention treats <2 as
    #     excellent, 5+ as concerning; 10 saturates anything truly bad).
    headers: Dict[str, dict] = {
        "mean_target_coverage": {
            "title": "Mean target cov.",
            "description": "Mean coverage depth across target bases (HQ non-dup)",
            "min": 0,
            "suffix": "X",
            "format": "{:,.1f}",
            "scale": "GnBu",
        },
        "frac_selected_bases": {
            "title": "% Selected",
            "description": "Fraction of aligned bases that are on or near a bait",
            "min": 50,
            "max": 100,
            "suffix": "%",
            "format": "{:,.1f}",
            "scale": "RdYlGn",
            "modify": lambda x: x * 100.0,
        },
        "frac_uncovered_targets": {
            "title": "% Uncovered targets",
            "description": "Fraction of raw (pre-merge) targets with zero coverage",
            "min": 0,
            "max": 10,
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "OrRd",
            "modify": lambda x: x * 100.0,
        },
        "hs_penalty_20x": {
            "title": "HS penalty 20X",
            "description": "Fold sequencing needed to reach 80% of targets at 20X",
            "min": 1,
            "max": 12,
            "format": "{:,.2f}",
            "scale": "OrRd",
            "hidden": True,
        },
        "at_dropout": {
            "title": "AT dropout",
            "description": "Under-representation of AT-rich targets (lower is better)",
            "min": 0,
            "max": 10,
            "format": "{:,.2f}",
            "scale": "OrRd",
            "hidden": True,
        },
        "gc_dropout": {
            "title": "GC dropout",
            "description": "Under-representation of GC-rich targets (lower is better)",
            "min": 0,
            "max": 10,
            "format": "{:,.2f}",
            "scale": "OrRd",
            "hidden": True,
        },
    }
    for cov in gs_covs:
        headers[f"frac_target_bases_{cov}x"] = {
            "title": f"Targets ≥ {cov}X",
            "description": f"Fraction of target bases with coverage ≥ {cov}X",
            "min": 80,
            "max": 100,
            "suffix": "%",
            "format": "{:,.1f}",
            "scale": "RdYlGn",
            "modify": lambda x: x * 100.0,
        }
    module.general_stats_addcols(data_by_sample, headers, namespace="hybcap")

    # Full metrics table — most columns hidden by default; the user can toggle them.
    table_data: Dict[str, Dict[str, float]] = {}
    for s_name, row in data_by_sample.items():
        merged = dict(row)
        merged["panel_name"] = panel_by_sample.get(s_name, "")
        table_data[s_name] = merged

    # Helper: percent column from a 0..1 fraction.
    # `direction` controls the colour gradient:
    #   "higher_better" -> RdYlGn (red at low end, green at high)
    #   "lower_better"  -> OrRd   (white-ish at 0, red at high)
    #   "neutral"       -> GnBu   (single hue, magnitude only)
    # `gmin`/`gmax` are the *gradient* bounds in percent units (after modify).
    def _pct(
        title: str,
        description: str,
        *,
        direction: str = "higher_better",
        hidden: bool = True,
        decimals: int = 2,
        gmin: float = 0,
        gmax: float = 100,
    ) -> dict:
        scale = {"higher_better": "RdYlGn", "lower_better": "OrRd", "neutral": "GnBu"}[direction]
        return {
            "title": title,
            "description": description,
            "min": gmin,
            "max": gmax,
            "suffix": "%",
            "format": f"{{:,.{decimals}f}}",
            "modify": lambda x: x * 100.0,
            "scale": scale,
            "hidden": hidden,
        }

    # An hs_penalty cell is "1.0 = perfect, higher = worse". Cap gradient at 12;
    # values up to that are routinely seen at deeper-coverage thresholds.
    def _hs_penalty(title: str, description: str, hidden: bool = True) -> dict:
        return {
            "title": title,
            "description": description,
            "min": 1,
            "max": 12,
            "format": "{:,.2f}",
            "scale": "OrRd",
            "hidden": hidden,
        }

    table_headers: Dict[str, dict] = {
        "panel_name": {"title": "Panel", "description": "Bait/target panel name", "hidden": True},
        "total_reads": {
            "title": "Total reads",
            "format": "{:,.0f}",
            "min": 0,
            "scale": "GnBu",
        },
        "frac_selected_bases": _pct(
            "% Selected",
            "Fraction of aligned bases that are on or near a bait",
            direction="higher_better",
            hidden=False,
            gmin=50,
        ),
        "mean_bait_coverage": {
            "title": "Mean bait cov.",
            "suffix": "X",
            "format": "{:,.1f}",
            "scale": "GnBu",
            "hidden": True,
        },
        "mean_target_coverage": {
            "title": "Mean target cov.",
            "suffix": "X",
            "format": "{:,.1f}",
            "scale": "GnBu",
        },
        "frac_uncovered_targets": _pct(
            "% Uncovered targets",
            "Fraction of raw (pre-merge) targets with zero coverage",
            direction="lower_better",
            hidden=False,
            gmax=10,
        ),
        "frac_exc_dupe": _pct(
            "% Exc dup", "Bases excluded as duplicates", direction="lower_better", gmax=30
        ),
        "frac_exc_mapq": _pct(
            "% Exc mapq", "Bases excluded due to low mapping quality", direction="lower_better", gmax=30
        ),
        "frac_exc_overlap": _pct(
            "% Exc overlap",
            "Bases excluded due to read-pair overlap clipping",
            direction="lower_better",
            gmax=30,
        ),
        "frac_exc_baseq": _pct(
            "% Exc baseq", "Bases excluded due to low base quality", direction="lower_better", gmax=30
        ),
        "frac_exc_off_target": _pct(
            "% Exc off-target", "Bases excluded as off-target (HQ)", direction="lower_better", gmax=60
        ),
        "frac_target_bases_1x": _pct(
            "% Targets ≥ 1X",
            "Fraction of target bases with coverage ≥ 1X",
            direction="higher_better",
            gmin=80,
        ),
        "frac_target_bases_20x": _pct(
            "% Targets ≥ 20X",
            "Fraction of target bases with coverage ≥ 20X",
            direction="higher_better",
            gmin=80,
        ),
        "frac_target_bases_50x": _pct(
            "% Targets ≥ 50X",
            "Fraction of target bases with coverage ≥ 50X",
            direction="higher_better",
            hidden=False,
            gmin=50,
        ),
        "frac_target_bases_100x": _pct(
            "% Targets ≥ 100X",
            "Fraction of target bases with coverage ≥ 100X",
            direction="higher_better",
            gmin=20,
        ),
        "frac_target_bases_250x": _pct(
            "% Targets ≥ 250X",
            "Fraction of target bases with coverage ≥ 250X",
            direction="higher_better",
        ),
        "frac_target_bases_500x": _pct(
            "% Targets ≥ 500X",
            "Fraction of target bases with coverage ≥ 500X",
            direction="higher_better",
        ),
        "at_dropout": {
            "title": "AT dropout",
            "description": "Under-representation of AT-rich targets (lower is better)",
            "min": 0,
            "max": 10,
            "format": "{:,.2f}",
            "scale": "OrRd",
            "hidden": True,
        },
        "gc_dropout": {
            "title": "GC dropout",
            "description": "Under-representation of GC-rich targets (lower is better)",
            "min": 0,
            "max": 10,
            "format": "{:,.2f}",
            "scale": "OrRd",
            "hidden": True,
        },
        "frac_usable_bases_on_target": _pct(
            "% Usable on target",
            "Fraction of total bases that are on-target (HQ, non-dup numerator)",
            direction="higher_better",
            gmin=20,
        ),
        "hs_library_size": {
            "title": "Library size (est.)",
            "description": "Estimated library size from selected pairs (Lander-Waterman)",
            "format": "{:,.0f}",
            "scale": "GnBu",
            "hidden": True,
        },
        "hs_penalty_20x": _hs_penalty(
            "HS penalty 20X", "Fold sequencing needed to reach 80% of targets at 20X"
        ),
        "hs_penalty_50x": _hs_penalty(
            "HS penalty 50X",
            "Fold sequencing needed to reach 80% of targets at 50X",
            hidden=False,
        ),
    }
    table_config = TableConfig(
        id=f"{module.anchor}_hybcap_table",
        title=f"{module.name}: Hybrid capture metrics",
    )
    module.add_section(
        name="Hybrid capture metrics",
        anchor=f"{module.anchor}-hybcap-table",
        description="Selected hybrid capture metrics. Many additional columns are available — use the column toggle to show/hide.",
        plot=table.plot(table_data, table_headers, table_config),
    )

    # Target coverage curve: % of target bases at >=Nx for the standard riker thresholds.
    curve_data: Dict[str, Dict[int, float]] = {}
    for s_name, row in data_by_sample.items():
        curve_data[s_name] = {}
        for cov in TARGET_COVERAGE_LEVELS:
            key = f"frac_target_bases_{cov}x"
            if key in row:
                curve_data[s_name][cov] = row[key] * 100.0

    if any(curve_data.values()):
        pconfig = {
            "id": f"{module.anchor}_hybcap_target_coverage",
            "title": f"{module.name}: Hybcap target coverage",
            "ylab": "% of target bases",
            "xlab": "Coverage depth",
            "ymin": 0,
            "ymax": 100,
            "x_decimals": False,
            "categories": True,
            "tt_label": "<b>≥{point.x}X</b>: {point.y:.1f}%",
        }
        module.add_section(
            name="Hybcap target coverage",
            anchor=f"{module.anchor}-hybcap-target-coverage",
            description="Cumulative fraction of target bases at the standard depth thresholds emitted by riker's `hybcap` tool.",
            plot=linegraph.plot(curve_data, pconfig),
        )

    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_hybcap")
    return set(data_by_sample.keys())
