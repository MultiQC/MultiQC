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
#
# CLIPPING: with `filterMode: "none"` (see `layout_option` below), dataZoom only moves the
# axis view range and never removes out-of-range items from `data`, so a bin straddling the
# zoomed edge (or, at full view, a bin near the grid boundary given the 1px pad above) still
# reaches renderItem with its untouched geometry. The SVG renderer used here does not
# auto-clip custom series to the grid the way it clips built-in series (bar/line/etc), so an
# unclipped rect paints straight over the y-axis sample labels and the plot title. Rather
# than depend on a global `echarts.graphic.clipRectByRect` (which the SSR path below cannot
# assume is in scope: `static_export.py` executes this body via MiniRacer against a fresh V8
# context, `params.coordSys` and `api.*` are always provided by the calling renderItem
# machinery itself), the rect is intersected against `params.coordSys` (the cartesian grid's
# pixel rect: `{x, y, width, height}`) by hand: shrink to the overlap, or draw nothing
# (`return;`, i.e. `undefined`, which ECharts treats as "no element for this data item") when
# the rect falls fully outside the grid.
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
    "var rx0 = p0[0];"
    "var ry0 = p0[1] - height / 2;"
    "var rx1 = rx0 + width;"
    "var ry1 = ry0 + height;"
    "var grid = params.coordSys;"
    "var gx0 = grid.x;"
    "var gy0 = grid.y;"
    "var gx1 = gx0 + grid.width;"
    "var gy1 = gy0 + grid.height;"
    "var cx0 = Math.max(rx0, gx0);"
    "var cy0 = Math.max(ry0, gy0);"
    "var cx1 = Math.min(rx1, gx1);"
    "var cy1 = Math.min(ry1, gy1);"
    "if (cx1 <= cx0 || cy1 <= cy0) return;"
    "return {"
    "type: 'rect',"
    "shape: { x: cx0, y: cy0, width: cx1 - cx0, height: cy1 - cy0 },"
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
    # Both axes' "inside" dataZoom rely on the series' `encode` (see series() below) to
    # filter by the right item dimensions (x: [start, end], y: row); without it, the y
    # dataZoom silently dropped any bin whose `end` value exceeded the sample count
    # (Item 1's regression).
    zoom = zoom_option(x=True, y=True)
    # BUG FIX: default dataZoom `filterMode` ("filter") REMOVES any data item not fully
    # within the zoomed window from the series before renderItem runs, so a bin only
    # partially covered by the zoom window (e.g. positions 50-100 zoomed to 30-60)
    # vanished entirely instead of showing its visible portion (Plotly has no such bug,
    # since its heatmap is a raster image, not per-bin marks). "none" makes dataZoom only
    # change the axis view range and leaves the series' `data` untouched, so renderItem
    # still draws every bin and the grid clips whatever falls outside the visible range
    # (same fix as violin.py/violin.js's dataZoom, for the same underlying ECharts
    # behavior). Scoped to THIS plot's dataZoom entries only, not `zoom_option`'s shared
    # default, since other plot types' series rely on the default "filter" behavior.
    # Both the "inside" dataZoom entries AND the toolbox box-select `dataZoom` feature
    # (which independently creates its own dataZoom instance on drag) need the override.
    for dz in zoom["dataZoom"]:
        dz["filterMode"] = "none"
    zoom["toolbox"]["feature"]["dataZoom"]["filterMode"] = "none"
    option.update(zoom)
    # Item 3: positions start at 1, xAxis.min above already fixes the axis's own scale
    # extent there, and an "inside" dataZoom's start/end window can never exceed the
    # axis's configured min/max, so pan/zoom can't show x < 1 (let alone < 0).

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
            # REGRESSION FIX (Item 1): without an explicit `encode`, ECharts defaults a
            # cartesian custom series' dimension-to-axis mapping to dim0 -> x, dim1 -> y
            # (positional inference), i.e. it would treat `end` (item index 1) as the
            # Y VALUE. That's invisible for renderItem itself (it reads api.value(0..6)
            # directly, unaffected by encode), but `zoom_option`'s "inside" dataZoom on
            # the yAxis uses this inferred encoding to decide which items are "in range"
            # (filterMode "filter"): any bin whose `end` exceeded the sample count (e.g.
            # end=49 on a report with 34 samples/category indices 0..33) was silently
            # dropped from the whole series, even with no zoom applied (the dataZoom's
            # default range is the full 0..100%, but it still filters by this bogus
            # per-item comparison). That's why long-read samples' wide bins (position
            # >= n_samples) vanished after `zoom_option(y=True)` was added: this encode
            # makes the y dataZoom filter by `row` (item index 2, the real sample axis
            # value) instead, and the x dataZoom filter by the [start, end] pair, so a
            # bin is kept whenever start OR end falls in the zoomed window (no gaps at
            # the window edges when a bin straddles the boundary).
            "encode": {"x": [0, 1], "y": 2},
            # Item 4: this heatmap's primary interaction is click-to-drill-down (see
            # afterPlotCreated()/openDrilldown() in the JS twin), so a pointer cursor is
            # the correct affordance over the rects, not ECharts' default zoom/crosshair
            # cursor. `cursor` is a top-level series option (drawn from every rect this
            # series renders), so this only affects the seqcontent plot's own series,
            # never other echarts plot types.
            "cursor": "pointer",
        }
    ]


def axis_data(dataset: Dataset, pconfig: SeqContentConfig) -> List[Tuple[str, List[str]]]:
    """`[("yAxis", dataset.samples)]`; the xAxis is a value axis, no category data."""
    return [("yAxis", dataset.samples)]


def mark_count(dataset: Dataset) -> int:
    """Total bins across all samples (not `n_samples * max_bp`, see the module docstring)."""
    return sum(len(bins) for bins in dataset.rows)
