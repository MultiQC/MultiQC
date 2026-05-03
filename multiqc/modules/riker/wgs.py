"""Parse riker `wgs` outputs (wgs-metrics.txt + wgs-coverage.txt)."""

import logging
from collections import defaultdict
from typing import Dict, List

from multiqc import config
from multiqc.plots import bargraph, linegraph, table
from multiqc.plots.bargraph import BarPlotConfig
from multiqc.plots.table import TableConfig

from .util import read_tsv, to_float, to_int

log = logging.getLogger(__name__)

# Coverage thresholds emitted by riker as `frac_bases_at_{N}x`.
WGS_COVERAGE_LEVELS = [1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 100]


def parse_reports(module):
    parsed_samples = _parse_metrics(module)
    parsed_samples |= _parse_coverage_histogram(module)
    return parsed_samples


def _parse_metrics(module) -> set:
    data_by_sample: Dict[str, Dict[str, float]] = {}

    for f in module.find_log_files("riker/wgs_metrics", filehandles=True):
        for row in read_tsv(f["f"]):
            sample = row.pop("sample", None)
            if not sample:
                continue
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="wgs_metrics")
            data_by_sample[s_name] = {col: to_float(val) for col, val in row.items()}

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    riker_config = getattr(config, "riker_config", {})
    user_covs = riker_config.get("general_stats_target_coverage", [])
    if isinstance(user_covs, list) and user_covs:
        gs_covs: List[int] = [int(c) for c in user_covs if int(c) in WGS_COVERAGE_LEVELS]
    else:
        gs_covs = [30]

    # Coloring rationale:
    #   - mean_coverage / median_coverage: experiment-dependent, single-hue GnBu.
    #   - frac_excluded_total: lower is better, gradient capped at 30% (above
    #     that the run is in trouble).
    #   - frac_bases_at_Nx: higher is better, gradient over 80-100% so a typical
    #     well-run sample doesn't paint the same shade as a marginal one.
    headers: Dict[str, dict] = {
        "median_coverage": {
            "title": "Median cov.",
            "description": "Median coverage depth",
            "min": 0,
            "suffix": "X",
            "format": "{:,.1f}",
            "scale": "GnBu",
        },
        "mean_coverage": {
            "title": "Mean cov.",
            "description": "Mean coverage depth",
            "min": 0,
            "suffix": "X",
            "format": "{:,.1f}",
            "scale": "GnBu",
            "hidden": True,
        },
        "frac_excluded_total": {
            "title": "% Excluded",
            "description": "Total fraction of bases excluded for any reason",
            "min": 0,
            "max": 30,
            "suffix": "%",
            "format": "{:,.1f}",
            "scale": "OrRd",
            "modify": lambda x: x * 100.0,
        },
    }
    for cov in gs_covs:
        headers[f"frac_bases_at_{cov}x"] = {
            "title": f"Bases ≥ {cov}X",
            "description": f"Fraction of genome territory with coverage ≥ {cov}X",
            "min": 80,
            "max": 100,
            "suffix": "%",
            "format": "{:,.1f}",
            "scale": "RdYlGn",
            "modify": lambda x: x * 100.0,
        }
    module.general_stats_addcols(data_by_sample, headers, namespace="wgs")

    # Detailed metrics table.
    # `direction` controls the colour gradient (see hybcap.py for the same helper):
    #   "higher_better" -> RdYlGn (gradient over gmin..gmax = 80..100 by default)
    #   "lower_better"  -> OrRd   (gradient over 0..gmax = 0..20 by default)
    #   "neutral"       -> GnBu (single hue, magnitude only)
    def _pct(
        title: str,
        description: str,
        *,
        direction: str = "higher_better",
        hidden: bool = True,
        gmin: float = 80,
        gmax: float = 100,
    ) -> dict:
        scale = {"higher_better": "RdYlGn", "lower_better": "OrRd", "neutral": "GnBu"}[direction]
        return {
            "title": title,
            "description": description,
            "min": gmin,
            "max": gmax,
            "suffix": "%",
            "format": "{:,.2f}",
            "modify": lambda x: x * 100.0,
            "scale": scale,
            "hidden": hidden,
        }

    table_headers: Dict[str, dict] = {
        "mean_coverage": {
            "title": "Mean cov.",
            "description": "Mean coverage depth",
            "suffix": "X",
            "format": "{:,.2f}",
            "scale": "GnBu",
            "min": 0,
        },
        "sd_coverage": {
            "title": "SD cov.",
            "description": "Standard deviation of coverage",
            "suffix": "X",
            "format": "{:,.2f}",
            "scale": "GnBu",
            "min": 0,
        },
        "median_coverage": {
            "title": "Median cov.",
            "description": "Median coverage depth",
            "suffix": "X",
            "format": "{:,.2f}",
            "scale": "GnBu",
            "min": 0,
        },
        "frac_excluded_mapq": _pct(
            "% Exc mapq", "Bases excluded due to low mapping quality", direction="lower_better", gmin=0, gmax=20
        ),
        "frac_excluded_dupe": _pct(
            "% Exc dupe", "Bases excluded as duplicates", direction="lower_better", gmin=0, gmax=30
        ),
        "frac_excluded_unpaired": _pct(
            "% Exc unpaired",
            "Bases excluded for unpaired / mate-unmapped reads",
            direction="lower_better",
            gmin=0,
            gmax=10,
        ),
        "frac_excluded_baseq": _pct(
            "% Exc baseq",
            "Bases excluded due to low base quality",
            direction="lower_better",
            gmin=0,
            gmax=10,
        ),
        "frac_excluded_overlap": _pct(
            "% Exc overlap",
            "Bases excluded due to read-pair overlap",
            direction="lower_better",
            gmin=0,
            gmax=10,
        ),
        "frac_excluded_capped": _pct(
            "% Exc capped",
            "Bases excluded for exceeding the coverage cap",
            direction="lower_better",
            gmin=0,
            gmax=5,
        ),
        "frac_excluded_total": _pct(
            "% Exc total",
            "Total fraction of bases excluded",
            direction="lower_better",
            hidden=False,
            gmin=0,
            gmax=30,
        ),
        "frac_bases_at_1x": _pct(
            "% Bases ≥ 1X",
            "Fraction of genome territory with coverage ≥ 1X",
            hidden=False,
            gmin=80,
        ),
        "frac_bases_at_5x": _pct(
            "% Bases ≥ 5X",
            "Fraction of genome territory with coverage ≥ 5X",
            hidden=False,
            gmin=80,
        ),
        "frac_bases_at_10x": _pct(
            "% Bases ≥ 10X",
            "Fraction of genome territory with coverage ≥ 10X",
            gmin=80,
        ),
        "frac_bases_at_20x": _pct(
            "% Bases ≥ 20X",
            "Fraction of genome territory with coverage ≥ 20X",
            gmin=80,
        ),
        "frac_bases_at_30x": _pct(
            "% Bases ≥ 30X",
            "Fraction of genome territory with coverage ≥ 30X",
            hidden=False,
            gmin=80,
        ),
    }
    table_config = TableConfig(
        id=f"{module.anchor}_wgs_table",
        title=f"{module.name}: WGS coverage metrics",
    )
    module.add_section(
        name="WGS coverage metrics",
        anchor=f"{module.anchor}-wgs-table",
        description="Per-sample WGS coverage and exclusion metrics from riker's `wgs` tool.",
        plot=table.plot(data_by_sample, table_headers, table_config),
    )

    # Bar plot: fraction of bases excluded for each reason.
    excl_keys = [
        ("frac_excluded_mapq", "Low mapping quality"),
        ("frac_excluded_dupe", "Duplicates"),
        ("frac_excluded_unpaired", "Unpaired / mate unmapped"),
        ("frac_excluded_baseq", "Low base quality"),
        ("frac_excluded_overlap", "Read-pair overlap"),
        ("frac_excluded_capped", "Coverage cap"),
    ]
    excl_data: Dict[str, Dict[str, float]] = {}
    for s_name, row in data_by_sample.items():
        excl_data[s_name] = {col: row.get(col, 0.0) * 100.0 for col, _ in excl_keys}
    excl_cats = {col: {"name": label} for col, label in excl_keys}
    excl_config = BarPlotConfig(
        id=f"{module.anchor}_wgs_excluded_bases",
        title=f"{module.name}: WGS excluded bases",
        ylab="% Bases excluded",
        ymax=100,
        cpswitch=False,
    )
    module.add_section(
        name="WGS excluded bases",
        anchor=f"{module.anchor}-wgs-excluded",
        description="Fraction of bases excluded from WGS coverage calculations, by exclusion reason.",
        plot=bargraph.plot(excl_data, excl_cats, excl_config),
    )

    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_wgs")
    return set(data_by_sample.keys())


def _parse_coverage_histogram(module) -> set:
    # data_by_sample[sample][depth] = bases (count at exactly that depth)
    data_by_sample: Dict[str, Dict[int, int]] = defaultdict(dict)

    for f in module.find_log_files("riker/wgs_coverage", filehandles=True):
        rows_by_sample: Dict[str, Dict[int, int]] = defaultdict(dict)
        for row in read_tsv(f["f"]):
            sample = row.get("sample")
            if not sample:
                continue
            try:
                depth = to_int(row["depth"])
                bases = to_int(row.get("bases", "0"))
            except (KeyError, ValueError):
                continue
            rows_by_sample[sample][depth] = bases

        for sample, by_depth in rows_by_sample.items():
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="wgs_coverage")
            data_by_sample[s_name] = by_depth

    data_by_sample = dict(module.ignore_samples(data_by_sample))
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    # Cut the histogram tail at 99% of bases to avoid extremely long tails dominating the plot.
    riker_config = getattr(config, "riker_config", {})
    max_cov = riker_config.get("wgs_histogram_max_cov")
    if max_cov is None:
        max_cov = 10
        for hist in data_by_sample.values():
            total = float(sum(hist.values()))
            if total <= 0:
                continue
            running = 0.0
            for depth, count in sorted(hist.items()):
                running += count
                if running > total * 0.99:
                    max_cov = max(depth, max_cov)
                    break

    trimmed: Dict[str, Dict[int, int]] = {}
    for s_name, hist in data_by_sample.items():
        trimmed[s_name] = {d: c for d, c in hist.items() if d <= max_cov}

    pconfig = {
        "id": f"{module.anchor}_wgs_coverage_histogram",
        "title": f"{module.name}: WGS coverage",
        "ylab": "Number of bases",
        "xlab": "Coverage depth",
        "x_decimals": False,
        "ymin": 0,
        "smooth_points": 500,
        "smooth_points_sumcounts": True,
        "tt_label": "<b>{point.x}X</b>: {point.y:,.0f} bases",
    }
    module.add_section(
        name="WGS coverage",
        anchor=f"{module.anchor}-wgs-coverage",
        description="Per-base coverage histogram from riker's `wgs` tool. The tail beyond 99% of bases is hidden.",
        plot=linegraph.plot(trimmed, pconfig),
    )

    return set(data_by_sample.keys())
