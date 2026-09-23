"""Family size histograms: `group --family-size-histogram`, `group --metrics` (`<prefix>.family_sizes.txt`)
and `dedup --family-size-histogram`. The columns are identical to fgbio GroupReadsByUmi's histogram, so the
fgumi search pattern is ordered to claim these files ahead of the fgbio module."""

import logging
from typing import Dict, Optional, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

from .schemas import FamilySizeMetric
from .util import drop_none, load_rows, pct, sample_name

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    counts: Dict[str, Dict[int, float]] = {}
    percent: Dict[str, Dict[int, Optional[float]]] = {}
    cumulative: Dict[str, Dict[int, Optional[float]]] = {}
    for f in module.find_log_files("fgumi/family_sizes"):
        rows = load_rows(f, FamilySizeMetric)
        if rows is None:
            continue
        s_name = sample_name(module, f)
        if s_name in counts:
            log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
        module.add_data_source(f, s_name)
        module.add_software_version(None, s_name)
        counts[s_name] = {r.family_size: r.count for r in rows}
        percent[s_name] = {r.family_size: pct(r.fraction) for r in rows}
        cumulative[s_name] = {r.family_size: pct(r.fraction_gt_or_eq_family_size) for r in rows}

    counts = module.ignore_samples(counts)
    if not counts:
        return set()
    percent = {s: percent[s] for s in counts}
    cumulative = {s: cumulative[s] for s in counts}

    module.add_section(
        name="Family sizes",
        anchor="fgumi-family-sizes",
        description="Distribution of UMI family sizes: the number of templates grouped into each molecule.",
        helptext="Written by `fgumi group` (`--family-size-histogram` or `--metrics`) and "
        "`fgumi dedup --family-size-histogram`. The cumulative tab shows the percentage of families of at least "
        "each size.",
        plot=linegraph.plot(
            [
                counts,
                {s: drop_none(v) for s, v in percent.items()},
                {s: drop_none(v) for s, v in cumulative.items()},
            ],
            {
                "id": "fgumi_family_sizes",
                "title": "fgumi: Family sizes",
                "xlab": "Family size (templates)",
                "xlog": True,
                "data_labels": [
                    {"name": "Counts", "ylab": "Families"},
                    {"name": "Percentages", "ylab": "% of families", "ysuffix": "%"},
                    {"name": "Cumulative", "ylab": "% of families at or above size", "ysuffix": "%"},
                ],
            },
        ),
    )
    module.write_data_file(counts, "multiqc_fgumi_family_sizes")
    return set(counts)
