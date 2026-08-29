"""
ECharts option builder for `PlotType.BAR` (`multiqc/plots/bargraph.py`).

MultiQC bar plots are horizontal: samples on the (inverted) category axis (y),
category values on the value axis (x). See `templates/default/src/js/plots/bar.js`
for the Plotly-JS equivalent this mirrors.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Tuple

from multiqc.plots.bargraph import BarPlotConfig, Dataset
from multiqc.plots.echarts.converter import convert_layout

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot


def layout_option(plot: "Plot[Any, Any]", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one bar-plot dataset, minus `series` and axis
    `data` arrays (samples are toolbox-dependent; the JS renderer fills `yAxis.data`).
    """
    option = convert_layout(plot.layout, dataset.layout)

    option["yAxis"]["type"] = "category"
    option["yAxis"]["inverse"] = True
    option["tooltip"]["trigger"] = "axis"

    # JS reads this to decide whether series stack (default) or sit side by side.
    option["_mqc"] = {"barmode": plot.layout.barmode}

    return option


def series(dataset: Dataset, pconfig: BarPlotConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    One `{"type": "bar"}` series per category. Mirrors the stacking logic in
    `BarPlot.initialize` (`multiqc/plots/bargraph.py`).
    """
    barmode = pconfig.stacking
    if barmode is None:  # legacy: non-default None means "group"
        barmode = "group"
    if barmode == "normal":  # legacy alias
        barmode = "relative"
    is_group = barmode == "group"

    result: List[Dict[str, Any]] = []
    for cat in dataset.cats:
        data = cat.data_pct if is_pct else cat.data
        s: Dict[str, Any] = {
            "type": "bar",
            "name": cat.name,
            "data": list(data),
            "itemStyle": {"color": cat.color},
            "barCategoryGap": "30%",
        }
        if not is_group:
            s["stack"] = "total"
        result.append(s)
    return result


def axis_data(dataset: Dataset) -> Tuple[str, List[str]]:
    """`("yAxis", sample names)`: bar plots are always horizontal, samples on `yAxis`."""
    return "yAxis", list(dataset.samples)


def mark_count(dataset: Dataset) -> int:
    return len(dataset.samples) * len(dataset.cats)
