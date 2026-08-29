"""
ECharts option builder for `PlotType.LINE` (`multiqc/plots/linegraph.py`).

See `templates/default/src/js/plots/line.js` for the Plotly-JS equivalent this mirrors,
and `templates/echarts/src/js/plots/line.js` for the interactive (toolbox-aware) JS
counterpart of `series()` below.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from multiqc.plots.echarts.converter import (
    bands_and_lines,
    convert_layout,
    echarts_dash,
    trailing_bands_lines_series,
    zoom_option,
)
from multiqc.plots.linegraph import Dataset, LinePlotConfig

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot


def _is_categorical(dataset: Dataset) -> bool:
    # Set by `linegraph.Dataset.create` when `pconfig.categories` (or a per-dataset
    # `dconfig["categories"]` override) is truthy; reused here instead of re-deriving it.
    return dataset.layout.get("xaxis", {}).get("type") == "category"


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
    categorical = _is_categorical(dataset)
    # yAxis is always the value axis; xAxis only is too when not categorical. Scale
    # (Plotly-style data-fitted autorange, not forced-0) doesn't apply to a category
    # axis, see `converter._convert_axis`.
    option = convert_layout(plot.layout, dataset.layout, scale_x=not categorical, scale_y=True)

    option["tooltip"]["trigger"] = "axis"
    # Plotly-style click+drag box-zoom on both axes (POLISH.md #17): both x and y carry
    # real data (a categorical x here is still a meaningful position, e.g. read-length
    # bins, not a sample list), see `zoom_option`.
    option.update(zoom_option(x=True, y=True))

    if categorical:
        option["xAxis"]["data"] = _categories(dataset)

    bands_lines = bands_and_lines(plot.pconfig)
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
    `trailing_bands_lines_series`. This is the SSR/get_option (non-toolbox) path; the
    interactive path is `EchartsLinePlot.buildSeries()` (`templates/echarts/src/js/plots/line.js`).

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
                "lineStyle": {"width": line.width, "type": echarts_dash(line.dash)},
                "itemStyle": {"color": line.color},
                "color": line.color,
            }
        )

    trailing = trailing_bands_lines_series(pconfig)
    if trailing:
        result.append(trailing)

    return result


def axis_data(dataset: Dataset, pconfig: LinePlotConfig) -> Optional[List[Tuple[str, List[str]]]]:
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
