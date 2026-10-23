"""Yield by input read pairs from `simplex-metrics` / `duplex-metrics` (fgbio CollectDuplexSeqMetrics plots #4a/#4b)."""

from typing import Dict, Optional, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

from .schemas import DuplexYieldMetric, SimplexYieldMetric
from .util import drop_none, finite, flatten, header_columns, load_rows, pct, register

_DUPLEX_TABS = [
    ("duplexes", "Duplexes (actual)", "Duplexes"),
    ("ideal_duplexes", "Duplexes (ideal)", "Duplexes"),
    ("ds_families", "Double-strand families", "Families"),
    ("ratio", "Actual / ideal duplexes", "Ratio"),
]
_SIMPLEX_TABS = [
    ("consensus_families", "Single-strand consensus families", "Families"),
    ("ss_families", "Single-strand families", "Families"),
    ("mean_family_size", "Mean single-strand family size", "Read pairs"),
]


def _ratio(actual: float, ideal: float) -> Optional[float]:
    """actual / ideal, or ``None`` when ideal is zero or non-finite (never a misleading 0)."""
    finite_ideal = finite(ideal)
    return finite(actual / finite_ideal) if finite_ideal else None


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    duplex: Dict[str, Dict[str, Dict[int, Optional[float]]]] = {}
    simplex: Dict[str, Dict[str, Dict[int, Optional[float]]]] = {}
    general: Dict[str, Dict[str, Optional[float]]] = {}
    for f in module.find_log_files("fgumi/yield"):
        if f["f"] is None:
            continue
        if DuplexYieldMetric.matches(header_columns(f)):
            rows = load_rows(f, DuplexYieldMetric)
            if not rows:
                continue
            s_name = register(module, f, section="duplex_yield")
            duplex[s_name] = {
                "duplexes": {r.read_pairs: r.ds_duplexes for r in rows},
                "ideal_duplexes": {r.read_pairs: finite(r.ds_families * r.ds_fraction_duplexes_ideal) for r in rows},
                "ds_families": {r.read_pairs: r.ds_families for r in rows},
                "ratio": {r.read_pairs: _ratio(r.ds_fraction_duplexes, r.ds_fraction_duplexes_ideal) for r in rows},
            }
            full = max(rows, key=lambda r: r.fraction)
            general.setdefault(s_name, {}).update(
                {"ds_duplexes": full.ds_duplexes, "ds_fraction_duplexes": pct(full.ds_fraction_duplexes)}
            )
        else:
            simplex_rows = load_rows(f, SimplexYieldMetric)
            if not simplex_rows:
                continue
            s_name = register(module, f, section="simplex_yield")
            simplex[s_name] = {
                "consensus_families": {r.read_pairs: r.ss_consensus_families for r in simplex_rows},
                "ss_families": {r.read_pairs: r.ss_families for r in simplex_rows},
                "mean_family_size": {r.read_pairs: finite(r.mean_ss_family_size) for r in simplex_rows},
            }
            simplex_full = max(simplex_rows, key=lambda r: r.fraction)
            general.setdefault(s_name, {})["ss_consensus_families"] = simplex_full.ss_consensus_families

    parsed: Set[str] = set()
    for kind, data, tabs in [("duplex", duplex, _DUPLEX_TABS), ("simplex", simplex, _SIMPLEX_TABS)]:
        data = module.ignore_samples(data)
        if not data:
            continue
        module.add_section(
            name=f"{kind.capitalize()} yield",
            anchor=f"fgumi-{kind}-yield",
            description=f"How {kind} yield grows with sequencing depth, from `fgumi {kind}-metrics`"
            + (" (or fgbio CollectDuplexSeqMetrics)" if kind == "duplex" else "")
            + " downsampling.",
            helptext="The input reads are downsampled to several fractions and the metrics recomputed at each, so the curves "
            "show how yield would grow with more sequencing. A curve that flattens means the library is saturating.",
            plot=linegraph.plot(
                [{s: drop_none(data[s][key]) for s in data} for key, _, _ in tabs],
                {
                    "id": f"fgumi_{kind}_yield",
                    "title": f"fgumi: {kind.capitalize()} yield",
                    "xlab": "Read pairs",
                    "data_labels": [{"name": name, "ylab": ylab} for _, name, ylab in tabs],
                },
            ),
        )
        module.write_data_file(flatten(data), f"multiqc_fgumi_{kind}_yield")
        parsed |= set(data)

    general = module.ignore_samples(general)
    if general:
        module.general_stats_addcols(
            {s: drop_none(v) for s, v in general.items()},
            {
                "ds_duplexes": {
                    "title": "Duplexes",
                    "scale": "YlGn",
                    "description": "Duplex consensus molecules",
                    "hidden": True,
                    "format": "{:,.0f}",
                },
                "ds_fraction_duplexes": {
                    "title": "% duplexes",
                    "scale": "RdYlGn",
                    "description": "Double-strand families that are duplexes",
                    "hidden": True,
                    "max": 100,
                    "suffix": "%",
                    "min": 0,
                },
                "ss_consensus_families": {
                    "title": "SS families",
                    "scale": "BuGn",
                    "description": "Single-strand consensus families",
                    "hidden": True,
                    "format": "{:,.0f}",
                },
            },
            namespace="yield",
        )
    return parsed
