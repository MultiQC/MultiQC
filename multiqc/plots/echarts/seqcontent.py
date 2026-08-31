"""
ECharts option builder for `PlotType.SEQCONTENT` (`multiqc/plots/seqcontent.py`).

Unlike the Plotly backend (a `go.Image` trace over a uniform per-bp pixel grid, see
`Dataset.create_figure`), ECharts draws ONE `custom` series `renderItem` rect per BIN,
with its true variable width (exactly what the old hand-written canvas drew): a run with
non-uniform bin widths across samples never needs an expanded `n_samples x max_bp` grid
here, which matters for long-read data. `series()` returns one `{"type": "custom"}` series
whose `data` items are `[start, end, row_idx, r, g, b, opacity]`, with r/g/b PRECOMPUTED in
Python (via `bin_rgb`, shared with the Plotly backend and the default template's JS, see
that function's own golden-contract docstring) so the SSR path never needs to run color
math in JS. `opacity` is the toolbox-highlight dim factor (always 1 here; only the live JS
side has highlight state).

`RENDER_ITEM_BODY` (the `__FN__` sentinel body executed by `static_export.py`'s walker,
verified to already support `series[].renderItem` with args `["params", "api"]`) is the
GOLDEN CONTRACT for its JS twin, the LIVE `renderItem` function built by
`EchartsSeqContentPlot.buildSeries()` in `templates/echarts/src/js/plots/seqcontent.js`
(Task 2.2): both must map a `[start, end, row, r, g, b, opacity]` item to the same rect
geometry, fill color, and opacity.
"""

from typing import Any, Dict, List, Tuple

from multiqc.plots.echarts.converter import convert_layout, zoom_option
from multiqc.plots.seqcontent import Dataset, SeqContentConfig, SeqContentPlot, bin_rgb

# GOLDEN CONTRACT: kept in lockstep with the live `renderItem` in
# `templates/echarts/src/js/plots/seqcontent.js` (Task 2.2). A bin covers columns
# `start..end` inclusive (1-based, matching `Dataset.create_figure`'s Plotly-side pixel
# fill), so its rect spans the data-space interval `[start, end + 1)` on the value xAxis;
# the yAxis is a category axis, one row per sample, so `row` is used as-is. The rect is
# padded 1px past its true data-space edges on both axes (width/height +1, height still
# centered on the row) so neighbouring bins overlap by ~1px instead of leaving a hairline
# seam between them, the way the old canvas renderer's solid per-pixel fill never showed
# gaps; no stroke is ever set, so there is no border to remove either. `opacity` (item
# index 6) is the toolbox-highlight dim factor computed by the live JS side's
# `buildSeries()`; the SSR/get_option path below always emits 1 (no highlight state
# exists server-side), but the field still has to exist so both twins share one rect
# shape.
RENDER_ITEM_BODY = (
    "var start = api.value(0);"
    "var end = api.value(1);"
    "var row = api.value(2);"
    "var r = api.value(3);"
    "var g = api.value(4);"
    "var b = api.value(5);"
    "var opacity = api.value(6);"
    "var p0 = api.coord([start, row]);"
    "var p1 = api.coord([end + 1, row]);"
    "var height = api.size([0, 1])[1] + 1;"
    "var width = p1[0] - p0[0] + 1;"
    "return {"
    "type: 'rect',"
    "shape: { x: p0[0], y: p0[1] - height / 2, width: width, height: height },"
    "style: { fill: 'rgb(' + r + ',' + g + ',' + b + ')', opacity: opacity }"
    "};"
)


def layout_option(plot: "SeqContentPlot", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one seqcontent dataset, minus `series` and yAxis
    `data` (the sample list, toolbox-dependent; filled in by `axis_data`/the JS renderer).
    """
    option = convert_layout(plot.layout, dataset.layout)

    option["xAxis"]["type"] = "value"
    option["xAxis"]["min"] = 1
    option["xAxis"]["max"] = dataset.max_bp + 1
    option["yAxis"]["type"] = "category"
    option["yAxis"]["inverse"] = True
    option["tooltip"]["trigger"] = "item"
    # No visualMap: colors are precomputed per-bin, not looked up from a value scale.
    # Unlike heatmap.py (a `heatmap`-type series, where the toolbox box-select brush was
    # confirmed empirically to never engage), this is a `custom` series with a real
    # numeric x-axis (bp position), so both axes get click+drag box-zoom like line/scatter.
    option.update(zoom_option(x=True, y=True))

    return option


def series(dataset: Dataset, pconfig: SeqContentConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    One `{"type": "custom"}` series, one rect per bin. This is the SSR/get_option
    (non-toolbox) path; the interactive path is `EchartsSeqContentPlot.buildSeries()`
    (`templates/echarts/src/js/plots/seqcontent.js`).

    `is_pct` is accepted for dispatch-signature parity with `bar.series`; seqcontent has
    no percentage switch, so it is unused here.

    Each data item is `[start, end, row_idx, r, g, b, opacity]`; `opacity` is always 1
    here since there is no toolbox highlight state server-side (see `RENDER_ITEM_BODY`'s
    docstring above).
    """
    data: List[List[int]] = []
    for row_idx, bins in enumerate(dataset.rows):
        for b in bins:
            r, g, bl = bin_rgb(b)
            data.append([b.start, b.end, row_idx, r, g, bl, 1])

    return [
        {
            "type": "custom",
            "coordinateSystem": "cartesian2d",
            "renderItem": {"__FN__": True, "body": RENDER_ITEM_BODY},
            "data": data,
        }
    ]


def axis_data(dataset: Dataset, pconfig: SeqContentConfig) -> List[Tuple[str, List[str]]]:
    """`[("yAxis", dataset.samples)]`; the xAxis is a value axis, no category data."""
    return [("yAxis", dataset.samples)]


def mark_count(dataset: Dataset) -> int:
    """Total bins across all samples (not `n_samples * max_bp`, see the module docstring)."""
    return sum(len(bins) for bins in dataset.rows)
