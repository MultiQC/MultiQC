"""`dedup --metrics` (one row per library plus `All Reads`) and `dedup --duplication-ladder`."""

import logging
from typing import Dict, Optional, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph

from .schemas import DeduplicationMetric, DuplicationLadderMetric
from .util import drop_none, flatten, iter_samples, load_rows, pct, register

ALL_READS = "All Reads"

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    return _metrics(module) | _ladder(module)


def _metrics(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, Optional[float]]] = {}
    for f in module.find_log_files("fgumi/dedup"):
        rows = load_rows(f, DeduplicationMetric)
        if not rows:
            continue
        # fgumi writes the aggregate row last; a real library could also be named "All Reads".
        total = next((r for r in reversed(rows) if r.library == ALL_READS), None)
        if total is None:
            log.warning(f"Skipping {f['fn']}: fgumi dedup metrics have no '{ALL_READS}' row")
            continue
        data[register(module, f)] = {
            "filtered_templates": total.filtered_templates,
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
        description="Unique, duplicate and filtered templates in `fgumi dedup` (all libraries combined).",
        helptext="Written by `fgumi dedup --metrics`, using the `All Reads` row. Filtered templates were dropped before "
        "duplicate marking (for example unmapped, low mapping quality, or a missing or N-containing UMI). % Dup in "
        "General Statistics is Picard's PERCENT_DUPLICATION.",
        plot=bargraph.plot(
            data,
            {
                "unique_templates": {"name": "Unique"},
                "duplicate_templates": {"name": "Duplicate"},
                "filtered_templates": {"name": "Filtered"},
            },
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
    module.write_data_file(flatten(data), "multiqc_fgumi_dedup")
    return set(data)


def _ladder(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, Dict[int, Optional[float]]]] = {}
    samples: Set[str] = set()
    # Multi-library rows are keyed "<sample> (<library>)"; iter_samples applies sample filters to the bare name.
    for s_name, rows in iter_samples(module, "fgumi/dedup_ladder", DuplicationLadderMetric):
        samples.add(s_name)
        libraries = sorted({r.library for r in rows})
        for library in libraries:
            key = s_name if len(libraries) == 1 else f"{s_name} ({library})"
            lib_rows = [r for r in rows if r.library == library]
            data[key] = {
                "cumulative": {r.templates_seen: pct(r.duplicate_fraction) for r in lib_rows},
                "window": {r.templates_seen: pct(r.window_duplicate_fraction) for r in lib_rows},
            }
    if not data:
        return set()
    module.add_section(
        name="Duplication ladder",
        anchor="fgumi-dedup-ladder",
        description="Cumulative and per-window duplicate rate as templates are processed in coordinate order.",
        helptext="Written by `fgumi dedup --duplication-ladder`, one line per library. Templates are counted in "
        "coordinate order, not in random order, so this is not a saturation curve: the rate at a point is the "
        "duplicate rate of the genome covered so far, and changes along the curve reflect regions with different "
        "duplicate rates.",
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
    module.write_data_file(flatten(data), "multiqc_fgumi_dedup_ladder")
    return samples
