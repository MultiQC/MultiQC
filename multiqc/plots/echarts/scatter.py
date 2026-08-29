"""
ECharts option builder for `PlotType.SCATTER` (`multiqc/plots/scatter.py`).

Series strategy: ONE `{"type": "scatter"}` series per dataset, with per-item styling
(`itemStyle.color`, `symbolSize`), rather than one series per point or per color-group.
The default template's Plotly `buildTraces()` (`templates/default/src/js/plots/scatter.js`)
builds one Plotly trace per point so every point can carry its own color/size/legend
entry independently; ECharts achieves the same per-point independence with per-item
`data` entries on a single series, which also keeps `mark_count` (point count) the
metric that decides the canvas-vs-svg renderer meaningful: thousands of series would
defeat that threshold's purpose, one series with thousands of data items does not.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from multiqc.plots.echarts.converter import convert_layout
from multiqc.plots.scatter import Dataset, PointT, ScatterConfig
from multiqc.utils.mqc_colour import color_to_rgb_string

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot

_DEFAULT_COLOR = "rgba(124, 181, 236, .5)"
_DEFAULT_SIZE = 10


def layout_option(plot: "Plot[Any, Any]", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one scatter-plot dataset, minus `series` (a single
    series built by the JS renderer from the raw points + toolbox state).
    """
    option = convert_layout(plot.layout, dataset.layout)

    option["tooltip"]["trigger"] = "item"
    option["dataZoom"] = [{"type": "inside"}, {"type": "slider"}]

    return option


def _point_item(point: PointT, default_color: str, default_size: Any) -> Dict[str, Any]:
    item: Dict[str, Any] = {
        "value": [point["x"], point["y"]],
        "name": point["name"],
        "symbolSize": point.get("marker_size", default_size),
        "itemStyle": {"color": color_to_rgb_string(str(point.get("color", default_color)))},
    }
    if "opacity" in point:
        item["itemStyle"]["opacity"] = point["opacity"]
    if "marker_line_width" in point:
        item["itemStyle"]["borderWidth"] = point["marker_line_width"]
    return item


def series(dataset: Dataset, pconfig: ScatterConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    One `{"type": "scatter"}` series holding every point in `dataset.points`, each with
    its own per-item style. This is the SSR/get_option (non-toolbox) path; the
    interactive path is `EchartsScatterPlot.buildSeries()`
    (`templates/echarts/src/js/plots/scatter.js`).

    `is_pct` is accepted for dispatch-signature parity with `bar.series`/`line.series`;
    scatter plots never enable the percentage switch (`Plot.initialize` only sets
    `add_pct_tab` for `PlotType.BAR`), so it is unused here.
    """
    marker = dataset.trace_params.get("marker", {})
    default_color = marker.get("color", _DEFAULT_COLOR)
    default_size = marker.get("size", _DEFAULT_SIZE)

    data = [_point_item(point, default_color, default_size) for point in dataset.points]

    return [
        {
            "type": "scatter",
            "name": dataset.label,
            "data": data,
        }
    ]


def axis_data(dataset: Dataset) -> Optional[Tuple[str, List[str]]]:
    """
    Always `None`: scatter plots use value (continuous) x/y axes, not a sample-name
    category axis, so there is no axis `data` array for the toolbox to fill in.
    """
    return None


def mark_count(dataset: Dataset) -> int:
    """One mark per point; this is the type most likely to trip the canvas threshold."""
    return len(dataset.points)
