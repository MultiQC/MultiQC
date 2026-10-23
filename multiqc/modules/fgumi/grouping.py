"""`group --grouping-metrics` / `<prefix>.grouping_metrics.txt` and `<prefix>.position_group_sizes.txt`."""

from typing import Dict, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph

from .schemas import PositionGroupSizeMetric, UmiGroupingMetric
from .util import drop_none, iter_samples, pct

_CATEGORIES = {
    "accepted_sam_records": {"name": "Accepted"},
    "discarded_non_pf": {"name": "Discarded: not passing filter"},
    "discarded_poor_alignment": {"name": "Discarded: poor alignment"},
    "discarded_ns_in_umi": {"name": "Discarded: Ns in UMI"},
    "discarded_umis_to_short": {"name": "Discarded: UMI too short"},
}


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    grouping: Dict[str, Dict[str, int]] = {
        s_name: rows[0].model_dump()
        for s_name, rows in iter_samples(module, "fgumi/grouping_metrics", UmiGroupingMetric)
        if rows
    }

    positions: Dict[str, Dict[int, float]] = {}
    positions_cumulative: Dict[str, Dict[int, float]] = {}
    for s_name, position_rows in iter_samples(module, "fgumi/position_group_sizes", PositionGroupSizeMetric):
        positions[s_name] = {r.position_group_size: r.count for r in position_rows}
        positions_cumulative[s_name] = drop_none(
            {r.position_group_size: pct(r.fraction_gt_or_eq_position_group_size) for r in position_rows}
        )

    if grouping:
        module.add_section(
            name="Grouping: primary records",
            anchor="fgumi-grouping-metrics",
            description="Primary alignment records accepted for grouping, and why the rest were discarded.",
            helptext="Written by `fgumi group` (`--grouping-metrics` or `--metrics`) or fgbio GroupReadsByUmi. "
            "Counts are primary records, "
            "matching fgbio's UmiGroupingMetric.",
            plot=bargraph.plot(
                grouping,
                _CATEGORIES,
                {"id": "fgumi_grouping_metrics", "title": "fgumi: Grouping", "ylab": "Primary records"},
            ),
        )
        module.write_data_file(grouping, "multiqc_fgumi_grouping_metrics")
    if positions:
        module.add_section(
            name="Grouping: position group sizes",
            anchor="fgumi-position-group-sizes",
            description="How many UMI families were found at each genomic position.",
            helptext="Written by `fgumi group --metrics` as `<prefix>.position_group_sizes.txt`.",
            plot=linegraph.plot(
                [positions, positions_cumulative],
                {
                    "id": "fgumi_position_group_sizes",
                    "title": "fgumi: Position group sizes",
                    "xlab": "Families at a position",
                    "xlog": True,
                    "data_labels": [
                        {"name": "Counts", "ylab": "Positions"},
                        {"name": "Cumulative", "ylab": "% of positions at or above size", "ysuffix": "%"},
                    ],
                },
            ),
        )
        module.write_data_file(positions, "multiqc_fgumi_position_group_sizes")
    return set(grouping) | set(positions)
