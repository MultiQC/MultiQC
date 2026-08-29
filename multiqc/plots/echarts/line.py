"""
ECharts option builder for `PlotType.LINE` (`multiqc/plots/linegraph.py`).

See `templates/default/src/js/plots/line.js` for the Plotly-JS equivalent this mirrors,
and `templates/echarts/src/js/plots/line.js` for the interactive (toolbox-aware) JS
counterpart of `series()` below.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from multiqc.plots.echarts.converter import convert_layout
from multiqc.plots.linegraph import Dataset, LinePlotConfig

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot

# `Series.dash` is already normalized to Plotly dash names by `convert_dash_style`
# (multiqc/plots/plot.py); ECharts `lineStyle.type` only knows "solid"/"dashed"/"dotted".
_DASH_MAP = {
    "solid": "solid",
    "dash": "dashed",
    "longdash": "dashed",
    "dot": "dotted",
    "dashdot": "dashed",
    "longdashdot": "dashed",
}


def _echarts_dash(dash: Optional[str]) -> str:
    if dash is None:
        return "solid"
    return _DASH_MAP.get(dash, "solid")


def _is_categorical(dataset: Dataset) -> bool:
    # Set by `linegraph.Dataset.create` when `pconfig.categories` (or a per-dataset
    # `dconfig["categories"]` override) is truthy; reused here instead of re-deriving it.
    return dataset.layout.get("xaxis", {}).get("type") == "category"


def _bands_and_lines(pconfig: LinePlotConfig) -> Dict[str, Any]:
    """
    Static markArea/markLine payload for `pconfig.x_bands`/`y_bands`/`x_lines`/`y_lines`
    (e.g. the fastqc per-base-quality green/yellow/red zones), mirroring the Plotly
    `shapes` built by `Plot._set_x_bands_and_range`/`_set_y_bands_and_range`
    (`multiqc/plots/plot.py`). These are static (independent of toolbox state), so this
    is computed once and reused by both `layout_option` (stashed in the skeleton under
    `_mqc.bandsLines` for the interactive JS path) and `series` (the SSR/get_option path).

    Per-dataset `data_labels` overrides of bands/lines are not supported yet (parity gap,
    same class of gap as bar's sample-groups gap); only the top-level `pconfig` fields
    are read.
    """
    mark_area_data: List[List[Dict[str, Any]]] = []
    for band in pconfig.y_bands or []:
        mark_area_data.append(
            [
                {"yAxis": band.from_, "itemStyle": {"color": band.color, "opacity": band.opacity}},
                {"yAxis": band.to},
            ]
        )
    for band in pconfig.x_bands or []:
        mark_area_data.append(
            [
                {"xAxis": band.from_, "itemStyle": {"color": band.color, "opacity": band.opacity}},
                {"xAxis": band.to},
            ]
        )

    mark_line_data: List[Dict[str, Any]] = []
    for line in pconfig.y_lines or []:
        mark_line_data.append(
            {
                "yAxis": line.value,
                "lineStyle": {"color": line.color, "width": line.width, "type": _echarts_dash(line.dash)},
                "label": {"formatter": line.label or ""},
            }
        )
    for line in pconfig.x_lines or []:
        mark_line_data.append(
            {
                "xAxis": line.value,
                "lineStyle": {"color": line.color, "width": line.width, "type": _echarts_dash(line.dash)},
                "label": {"formatter": line.label or ""},
            }
        )

    result: Dict[str, Any] = {}
    if mark_area_data:
        result["markArea"] = {"silent": True, "data": mark_area_data}
    if mark_line_data:
        result["markLine"] = {"silent": True, "symbol": "none", "data": mark_line_data}
    return result


def layout_option(plot: "Plot[Any, Any]", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one line-plot dataset, minus `series` (the JS
    renderer builds the toolbox-aware series list). `xAxis.type` ("category" vs "value")
    already comes out of `convert_layout` correctly, since `linegraph.Dataset.create`
    sets `dataset.layout["xaxis"]["type"]` up front for categorical plots.

    Deliberate deviation from the generic "no axis `data` arrays in the skeleton" rule:
    unlike bar's category axis (sample names, which toolbox hide/rename genuinely
    changes), a line plot's x-category labels (e.g. read-length bins) are shared,
    static values, identical for every sample and untouched by the toolbox. So they are
    safe, and simpler, to bake into the skeleton here rather than recompute in JS on
    every render; see `axis_data` below, which is a no-op for this type as a result.
    """
    option = convert_layout(plot.layout, dataset.layout)

    option["tooltip"]["trigger"] = "axis"
    option["dataZoom"] = [{"type": "inside"}, {"type": "slider"}]

    if _is_categorical(dataset):
        option["xAxis"]["data"] = _categories(dataset)

    bands_lines = _bands_and_lines(plot.pconfig)
    if bands_lines:
        option["_mqc"] = {"bandsLines": bands_lines}

    return option


def _categories(dataset: Dataset) -> List[str]:
    if not dataset.lines:
        return []
    return [str(x) for x, _ in dataset.lines[0].pairs]


def series(dataset: Dataset, pconfig: LinePlotConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    One `{"type": "line"}` series per `dataset.lines[i]`, plus (if configured) a trailing
    silent series carrying the static bands/lines markArea/markLine from
    `_bands_and_lines`. This is the SSR/get_option (non-toolbox) path; the interactive
    path is `EchartsLinePlot.buildSeries()` (`templates/echarts/src/js/plots/line.js`).

    `is_pct` is accepted for dispatch-signature parity with `bar.series`; line plots
    never enable the percentage switch (`Plot.initialize` only sets `add_pct_tab` for
    `PlotType.BAR`), so it is unused here.
    """
    categorical = _is_categorical(dataset)
    mode = dataset.trace_params.get("mode") or ""

    result: List[Dict[str, Any]] = []
    for line in dataset.lines:
        show_symbol = "markers" in mode or line.marker is not None
        data: List[Any] = [y for _, y in line.pairs] if categorical else [[x, y] for x, y in line.pairs]
        result.append(
            {
                "type": "line",
                "name": line.name,
                "data": data,
                "showSymbol": show_symbol,
                "smooth": False,
                "lineStyle": {"width": line.width, "type": _echarts_dash(line.dash)},
                "itemStyle": {"color": line.color},
                "color": line.color,
            }
        )

    bands_lines = _bands_and_lines(pconfig)
    if bands_lines:
        result.append(
            {
                "type": "line",
                "name": "",
                "data": [],
                "silent": True,
                "showSymbol": False,
                "tooltip": {"show": False},
                **bands_lines,
            }
        )

    return result


def axis_data(dataset: Dataset) -> Optional[Tuple[str, List[str]]]:
    """
    Always `None`: unlike bar's sample-name category axis, a line plot's x-category
    labels are static (see the comment in `layout_option`) and already baked into the
    skeleton there, so there is nothing left for `get_option`/JS to fill in here.
    """
    return None


def mark_count(dataset: Dataset) -> int:
    """One mark per line (path); plus one per point when markers are shown."""
    count = len(dataset.lines)
    mode = dataset.trace_params.get("mode") or ""
    if "markers" in mode:
        count += sum(len(line.pairs) for line in dataset.lines)
    return count
