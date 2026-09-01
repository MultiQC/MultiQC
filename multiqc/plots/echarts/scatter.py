"""
ECharts option builder for `PlotType.SCATTER` (`multiqc/plots/scatter.py`).

Series strategy: when no point carries a `group` field, ONE `{"type": "scatter"}` series
per dataset, with per-item styling (`itemStyle.color`, `symbolSize`), rather than one
series per point. The default template's Plotly `buildTraces()`
(`templates/default/src/js/plots/scatter.js`) builds one Plotly trace per point so every
point can carry its own color/size/legend entry independently; ECharts achieves the same
per-point independence with per-item `data` entries on a single series, which also keeps
`mark_count` (point count) the metric that decides the canvas-vs-svg renderer meaningful:
thousands of series would defeat that threshold's purpose, one series with thousands of
data items does not.

When points DO carry a `group` field (e.g. somalier's relatedness plot), that shortcut
would lose the per-group legend/toggle behaviour that is the entire point of grouping, so
`series()` instead builds ONE series PER GROUP (see `_grouped_series`), ordered per
`pconfig.groups`, mirroring Plotly's `legendgroup`/`inLegend` handling: ECharts gives each
series its own legend entry natively, and clicking it toggles that whole series. This is
still cheap because real grouped datasets have few groups (a handful of relatedness
categories, not thousands of points), so the canvas-threshold rationale above is
unaffected.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from multiqc.plots.echarts.converter import bands_and_lines, convert_layout, trailing_bands_lines_series, zoom_option
from multiqc.plots.scatter import Dataset, PointT, ScatterConfig
from multiqc.utils.mqc_colour import color_to_rgb_string

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot

_DEFAULT_COLOR = "rgba(124, 181, 236, .5)"
_DEFAULT_SIZE = 10

# Plotly marker symbol names (as used by `pconfig.marker_symbol`/`point["marker_symbol"]`)
# mapped to their closest ECharts built-in symbol. ECharts' built-in set is much smaller
# than Plotly's (no per-symbol open/dot variants, no exotic shapes), so this only covers
# the common names with a clean visual equivalent; anything else falls back to "circle",
# same as the interactive JS side (`scatter.js`'s `echartsSymbol()`).
_SYMBOL_MAP: Dict[str, str] = {
    "circle": "circle",
    "square": "rect",
    "diamond": "diamond",
    "triangle": "triangle",
    "triangle-up": "triangle",
}


def _echarts_symbol(marker_symbol: Optional[str]) -> str:
    if marker_symbol is None:
        return "circle"
    return _SYMBOL_MAP.get(marker_symbol, "circle")


def layout_option(plot: "Plot[Any, Any]", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one scatter-plot dataset, minus `series` (a single
    series built by the JS renderer from the raw points + toolbox state).
    """
    # Scatter has two meaningful value axes (no category axis), so both get Plotly-style
    # data-fitted autorange instead of ECharts' forced-0 default; see
    # `converter._convert_axis`.
    option = convert_layout(plot.layout_ir, dataset.layout, scale_x=True, scale_y=True)

    option["tooltip"]["trigger"] = "item"
    # Plotly-style click+drag box-zoom on both axes (POLISH.md #17): both axes carry
    # real numeric data, see `zoom_option`.
    option.update(zoom_option(x=True, y=True))

    # Scatter has two meaningful value axes (no category axis), so both dimensions
    # apply, same as line.
    bands_lines = bands_and_lines(plot.pconfig, dataset.dconfig)
    if bands_lines:
        option["_mqc"] = {"bandsLines": bands_lines}

    return option


def _point_item(point: PointT, default_color: str, default_size: Any) -> Dict[str, Any]:
    item: Dict[str, Any] = {
        "value": [point["x"], point["y"]],
        "name": point["name"],
        "symbolSize": point.get("marker_size", default_size),
        "symbol": _echarts_symbol(str(point["marker_symbol"]) if "marker_symbol" in point else None),
        "itemStyle": {"color": color_to_rgb_string(str(point.get("color", default_color)))},
    }
    if "opacity" in point:
        item["itemStyle"]["opacity"] = point["opacity"]
    if "marker_line_width" in point:
        item["itemStyle"]["borderWidth"] = point["marker_line_width"]
    return item


def _grouped_series(
    points: List[PointT], pconfig: ScatterConfig, default_color: str, default_size: Any
) -> List[Dict[str, Any]]:
    """
    One `{"type": "scatter"}` series per group, `name` = group name, so ECharts' legend
    shows one entry per group and toggling it hides/shows the whole group (ECharts' native
    per-series legend behaviour standing in for Plotly's `legendgroup`).

    Points are bucketed by `point["group"]`; a point with no group falls back to its own
    `point["name"]` as the bucket key, the same fallback the interactive JS
    (`scatter.js`'s `_buildGroupedSeries`) and the Plotly template's `buildTraces()`
    `inLegend` key use, so a point with no group still gets its own legend entry rather
    than being silently dropped from one.

    Buckets are ordered by first-seen point, then re-sorted by `pconfig.groups` when set
    (groups not listed there go last, keeping their relative first-seen order), mirroring
    the sort in the default template's Plotly `buildTraces()`.
    """
    order: List[str] = []
    buckets: Dict[str, List[PointT]] = {}
    for point in points:
        key = str(point.get("group") or point["name"])
        if key not in buckets:
            buckets[key] = []
            order.append(key)
        buckets[key].append(point)

    if pconfig.groups:
        group_order = pconfig.groups
        order.sort(key=lambda key: group_order.index(key) if key in group_order else len(group_order))

    return [
        {
            "type": "scatter",
            "name": key,
            "data": [_point_item(point, default_color, default_size) for point in buckets[key]],
        }
        for key in order
    ]


def series(dataset: Dataset, pconfig: ScatterConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    Scatter series for the SSR/get_option (non-toolbox) path; the interactive path is
    `EchartsScatterPlot.buildSeries()` (`templates/default/src/js/plots-echarts/scatter.js`).

    One `{"type": "scatter"}` series holding every point, each with its own per-item
    style, UNLESS any point carries a `group` field, in which case one series per group
    (see `_grouped_series`) so the legend and click-to-toggle mirror the Plotly template.

    `is_pct` is accepted for dispatch-signature parity with `bar.series`/`line.series`;
    scatter plots never enable the percentage switch (`Plot.initialize` only sets
    `add_pct_tab` for `PlotType.BAR`), so it is unused here.
    """
    marker = dataset.trace_params.get("marker", {})
    default_color = marker.get("color", _DEFAULT_COLOR)
    default_size = marker.get("size", _DEFAULT_SIZE)

    has_groups = any(point.get("group") for point in dataset.points)
    if has_groups:
        result: List[Dict[str, Any]] = _grouped_series(dataset.points, pconfig, default_color, default_size)
    else:
        result = [
            {
                "type": "scatter",
                "name": dataset.label,
                "data": [_point_item(point, default_color, default_size) for point in dataset.points],
            }
        ]

    trailing = trailing_bands_lines_series(pconfig, dataset.dconfig)
    if trailing:
        result.append(trailing)

    return result


def axis_data(dataset: Dataset, pconfig: ScatterConfig) -> Optional[List[Tuple[str, List[str]]]]:
    """
    Always `None`: scatter plots use value (continuous) x/y axes, not a sample-name
    category axis, so there is no axis `data` array for the toolbox to fill in.
    """
    return None


def mark_count(dataset: Dataset) -> int:
    """One mark per point; this is the type most likely to trip the canvas threshold."""
    return len(dataset.points)
