"""`dedup --metrics` (one row per library plus `All Reads`) and `dedup --duplication-ladder`."""

from typing import Dict, Optional, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph

from .schemas import DeduplicationMetric, DuplicationLadderMetric
from .util import drop_none, load_rows, pct, sample_name

ALL_READS = "All Reads"


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    return _metrics(module) | _ladder(module)


def _metrics(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, Optional[float]]] = {}
    for f in module.find_log_files("fgumi/dedup"):
        rows = load_rows(f, DeduplicationMetric)
        if not rows:
            continue
        total = next((r for r in rows if r.library == ALL_READS), rows[-1])
        s_name = sample_name(module, f)
        module.add_data_source(f, s_name)
        module.add_software_version(None, s_name)
        data[s_name] = {
            "unique_templates": total.unique_templates,
            "duplicate_templates": total.duplicate_templates,
            "percent_duplication": pct(total.percent_duplication),
            "estimated_library_size": total.estimated_library_size,
        }
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Deduplication",
        anchor="fgumi-dedup",
        description="Unique and duplicate templates found by `fgumi dedup` (all libraries combined).",
        plot=bargraph.plot(
            data,
            {"unique_templates": {"name": "Unique"}, "duplicate_templates": {"name": "Duplicate"}},
            {"id": "fgumi_dedup_templates", "title": "fgumi: Deduplication", "ylab": "Templates"},
        ),
    )
    module.general_stats_addcols(
        {
            s: drop_none(
                {"percent_duplication": v["percent_duplication"], "estimated_library_size": v["estimated_library_size"]}
            )
            for s, v in data.items()
        },
        {
            "percent_duplication": {
                "title": "% Dup",
                "description": "Duplication rate (Picard PERCENT_DUPLICATION)",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "scale": "OrRd",
            },
            "estimated_library_size": {
                "title": "Library size",
                "description": "Estimated library size",
                "format": "{:,.0f}",
                "hidden": True,
            },
        },
    )
    module.write_data_file(data, "multiqc_fgumi_dedup")
    return set(data)


def _ladder(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, Dict[int, Optional[float]]]] = {}
    for f in module.find_log_files("fgumi/dedup_ladder"):
        rows = load_rows(f, DuplicationLadderMetric)
        if rows is None:
            continue
        s_name = sample_name(module, f)
        module.add_data_source(f, s_name)
        module.add_software_version(None, s_name)
        libraries = sorted({r.library for r in rows})
        for library in libraries:
            key = s_name if len(libraries) == 1 else f"{s_name} ({library})"
            lib_rows = [r for r in rows if r.library == library]
            data[key] = {
                "cumulative": {r.templates_seen: pct(r.duplicate_fraction) for r in lib_rows},
                "window": {r.templates_seen: pct(r.window_duplicate_fraction) for r in lib_rows},
            }
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Duplication ladder",
        anchor="fgumi-dedup-ladder",
        description="Duplicate rate as templates accumulate: a flattening curve means the library is saturating.",
        helptext="Written by `fgumi dedup --duplication-ladder`.",
        plot=linegraph.plot(
            [
                {s: drop_none(v["cumulative"]) for s, v in data.items()},
                {s: drop_none(v["window"]) for s, v in data.items()},
            ],
            {
                "id": "fgumi_dedup_ladder",
                "title": "fgumi: Duplication ladder",
                "xlab": "Templates seen",
                "data_labels": [
                    {"name": "Cumulative", "ylab": "% duplicate", "ysuffix": "%"},
                    {"name": "Per window", "ylab": "% duplicate in window", "ysuffix": "%"},
                ],
            },
        ),
    )
    module.write_data_file(data, "multiqc_fgumi_dedup_ladder")
    return {key.split(" (", 1)[0] for key in data}
