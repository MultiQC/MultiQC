"""Family size histograms: `group --family-size-histogram`, `group --metrics` (`<prefix>.family_sizes.txt`)
and `dedup --family-size-histogram`. The columns are identical to fgbio GroupReadsByUmi's histogram, so both
modules find these files; `util.family_sizes_module` decides which one reports them."""

from typing import Dict, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

from .schemas import FamilySizeMetric
from .util import drop_none, family_sizes_evidence, family_sizes_module, found_by, iter_samples, pct


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    counts: Dict[str, Dict[int, float]] = {}
    percent: Dict[str, Dict[int, float]] = {}
    cumulative: Dict[str, Dict[int, float]] = {}
    evidence = family_sizes_evidence()
    for s_name, rows in iter_samples(
        module,
        "fgumi/family_sizes",
        FamilySizeMetric,
        skip=lambda f: found_by(f, "fgbio/groupreadsbyumi") and family_sizes_module(f, evidence) == "fgbio",
    ):
        counts[s_name] = {r.family_size: r.count for r in rows}
        percent[s_name] = drop_none({r.family_size: pct(r.fraction) for r in rows})
        cumulative[s_name] = drop_none({r.family_size: pct(r.fraction_gt_or_eq_family_size) for r in rows})
    if not counts:
        return set()

    module.add_section(
        name="Family sizes",
        anchor="fgumi-family-sizes",
        description="Distribution of UMI family sizes: the number of templates grouped into each molecule.",
        helptext="Written by `fgumi group` (`--family-size-histogram` or `--metrics`), fgbio GroupReadsByUmi, and "
        "`fgumi dedup --family-size-histogram`. The cumulative tab shows the percentage of families of at least "
        "each size.",
        plot=linegraph.plot(
            [counts, percent, cumulative],
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
