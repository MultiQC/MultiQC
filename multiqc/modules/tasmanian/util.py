"""Helpers shared by the tasmanian-diagnostics submodules"""

from typing import Dict, List

from multiqc import BaseMultiqcModule
from multiqc.plots import bargraph
from multiqc.types import SectionAlert


def add_category_bargraph(
    module: BaseMultiqcModule,
    counts: Dict[str, Dict[str, int]],
    *,
    name: str,
    anchor: str,
    description: str,
    helptext: str,
    ylab: str,
) -> None:
    """Add a stacked bar graph section of counts per category and sample, keeping the section if all counts are zero"""
    plot_data = {s_name: cats for s_name, cats in counts.items() if any(cats.values())}
    empty: List[str] = sorted(set(counts) - set(plot_data))
    alerts = (
        SectionAlert(
            message=f"**{len(empty)} sample{'s' if len(empty) != 1 else ''}** with no records hidden from this plot.",
            level="info" if plot_data else "warning",
            affected_samples=empty,
        )
        if empty
        else None
    )
    module.add_section(
        name=name,
        anchor=anchor,
        description=description,
        helptext=helptext,
        plot=bargraph.plot(
            plot_data,
            pconfig={
                "id": f"{anchor}-plot",
                "title": f"Tasmanian: {name}",
                "ylab": ylab,
                "cpswitch_counts_label": "Count",
            },
        )
        if plot_data
        else None,
        alerts=alerts,
    )
