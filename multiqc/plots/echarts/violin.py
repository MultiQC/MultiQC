"""
ECharts option builder for `PlotType.VIOLIN` (`multiqc/plots/violin.py`).

Ports `multiqc-echarts-exploration/scripts/06_violin_final.py` faithfully, then extends
it with PLOTLY-STYLE PER-ROW SUBPLOTS: one ECharts `grid` (with its own `xAxis`/`yAxis`
pair) per visible metric, stacked vertically, matching `multiqc/plots/violin.py`'s own
`Dataset.create_figure` (one Plotly subplot per metric, real x-axis values, compact
per-row height). Each row is drawn as a hand-built KDE polygon via a `custom` series
`renderItem`, plus a beeswarm `scatter` series with v6's native `jitter`, plus a `custom`
series per row for the min/max text annotations.

Per-row geometry: `_row_geometry` mirrors the same formula duplicated in
`templates/echarts/src/js/plots/violin.js` (`rowGeometry`), so a row's grid percentage
split is identical whether computed here (SSR/`get_option`) or in the browser
(interactive). Because grid percentages are resolved against whatever pixel height the
chart is actually given, this self-scales correctly for both the SSR/flat path (container
height comes from the shared Plotly-style `plot.layout.height`, unrelated to this file)
and the interactive path (the JS renderer sets a compact, row-count-scaled CSS height
itself; see `EchartsViolinPlot.buildSeries()`).

Since each row now gets its OWN real-valued x-axis, a row's KDE polygon/box/median/
scatter/annotation coordinates are drawn directly in REAL x (no more per-row 0..1
normalization); only the y coordinate stays a small symmetric offset around 0 (the
violin's density thickness), identical in scale to the old `ROW_HEIGHT`/`BOX_HALF_HEIGHT`/
`MEDIAN_HALF_HEIGHT` constants, just no longer shifted by a row index.

Two deliberate exceptions to the general "ECharts model->JSON contract"
(`multiqc-echarts-exploration/BUILD_PLAN.md`), both called out there by name:

1. VALUE-AXIS TRICK: each row's `yAxis` is a `type: "value"` axis (not `"category"`),
   spanning a fixed `[-0.5, 0.5]`, with a SINGLE tick at 0 whose label is the metric
   title, `verticalAlign: "middle"` (POLISH.md #1: the label is vertically centered on
   its row, since the tick sits at the exact vertical center of that row's own grid).
   `layout_option` therefore emits `yAxis.axisLabel.formatter` as a `__FN__` sentinel
   INSIDE the skeleton returned by `serialize()`, unlike every other type (whose
   skeletons never carry formatter functions, because the interactive JS renderer
   attaches its own). This is safe here because the row->title mapping never changes for
   hidden/renamed/highlighted samples (only for a change in which METRICS exist), and
   because the interactive JS path (`templates/echarts/src/js/plots/violin.js`) rebuilds
   `grid`/`xAxis`/`yAxis` wholesale from the LIVE metric list every render (which can
   differ from this module's SSR-time list, e.g. a table column toggled after page load)
   rather than patching this sentinel in place; only the SSR path (`static_export.py`'s
   `__FN__` walker) actually executes the sentinel this module emits.

2. Pre-computed polygons: `serialize()` additionally stores, per dataset, a `"violins"`
   key (`{metric: {"poly": [[x, y], ...], "range": [lo, hi]}}`) computed from ALL samples
   (not toolbox-filtered). KDE is too expensive to redo in JS on every render; the
   interactive renderer uses this polygon directly when no samples are hidden, and
   recomputes with the JS port of `kde()` only when samples are actually hidden. Because
   `poly`'s y no longer encodes a row offset (see above), the interactive renderer can
   reuse a precomputed polygon verbatim regardless of the row's live position, with no
   translation step needed (unlike the old design, where a metric's row shifting
   relative to Python's serialization order required shifting the cached polygon's y by
   hand). The inner Q1-Q3 box + median line (POLISH.md #10b) is NOT cached here: a sort +
   linear interpolation is cheap enough that the JS side (`buildSeries()`) just recomputes
   it from the currently-visible values on every render, same as the beeswarm points.

Trust boundary: every `{"__FN__": True, "body": "<js source>"}` sentinel emitted in this
module is built from typed, MultiQC-generated values (floats, sample names formatted with
`json.dumps`), never raw external input; `static_export.py` turns these bodies into real
JS functions via `new Function(...)`, which is only safe because of that provenance (see
the comment at `static_export._FN_WALKER_JS`).
"""

import json
import re
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Tuple

from multiqc.plots.echarts.box import _quantile
from multiqc.plots.echarts.converter import _clean_title_segment, _si_axis_formatter_body, convert_layout
from multiqc.plots.table_object import TableConfig
from multiqc.plots.violin import ColumnAnchor, Dataset, ViolinColumn
from multiqc.utils.mqc_colour import color_to_rgb_string

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot

ROW_HEIGHT = 0.42
N_BINS = 60
RANGE_PAD = 0.15

# Inner Q1-Q3 box + median line drawn over each violin (POLISH.md #10b), sized as
# fractions of ROW_HEIGHT so they scale with it; roughly matches Plotly's own default
# violin box width ratio (~0.25 of the violin's full height).
BOX_HALF_HEIGHT = 0.12
MEDIAN_HALF_HEIGHT = 0.16

# Fill opacity for the violin body (POLISH.md #10a). Matches Plotly's own dark-theme
# violin fill opacity (see templates/default/src/js/plots/violin.js's
# `isDarkMode ? 0.3 : 0.5`, applied client-side after Plotly renders); the SSR path here
# has no way to detect the viewer's theme, so both this module and the interactive JS
# path fix on the single dark-theme value for consistency (see module docstring).
FILL_ALPHA = 0.3
# Box fill is more solid than the violin body so the Q1-Q3 range reads clearly against it.
BOX_FILL_ALPHA = 0.55

# Default violin body color (POLISH.md #5): matches Plotly's own general-stats violin
# default, `fillcolor="#999999"` (see `multiqc/plots/violin.py`'s `Dataset.create`), kept
# as rgb(...) so _rgba() below applies uniformly. Only used when a metric has no
# configured `header.color`; a real per-metric color always takes priority (see
# `_violin_colors` below).
_DEFAULT_STROKE = "rgb(153,153,153)"
_DEFAULT_FILL = f"rgba(153,153,153,{FILL_ALPHA})"
_DEFAULT_BOX_FILL = f"rgba(153,153,153,{BOX_FILL_ALPHA})"
_SCATTER_COLOR = "rgba(30,50,80,0.85)"

# --- Per-row grid geometry: the cross-language contract --------------------------------
#
# GOLDEN CROSS-LANGUAGE CONTRACT: kept in lockstep with `rowGeometry()` in
# `templates/echarts/src/js/plots/violin.js` (same duplication pattern as `kde()`).
# `ROW_PX`/`EXTRA_PX`/`BOTTOM_PX` describe an IDEAL container height (`ROW_PX * n +
# EXTRA_PX + BOTTOM_PX`) that the interactive JS renderer actually builds (via a CSS
# height it sets itself, giving Plotly-comparable compact rows, POLISH.md #4); this
# module has no control over the real SSR/flat container height (that comes from the
# shared Plotly-style `plot.layout.height`, part of the data layer), but computing
# percentages against the SAME ideal target keeps row proportions identical between the
# two paths whenever the actual height happens to match, and merely rescales
# proportionally (not incorrectly) when it doesn't.
ROW_PX = 42.0
EXTRA_PX = 55.0  # space reserved for the chart title/subtitle above the first row
BOTTOM_PX = 10.0  # breathing room below the last row's x-axis tick labels
# Fraction of a row's slot handed to ECharts as the grid box (`containLabel: true` then
# shrinks the ACTUAL drawing area within that box to make room for the row's own x-axis
# tick labels, so this only needs to cover the tick labels plus a small gap to the next
# row, NOT also the violin height on top of that: 34px * (1 - 0.86) ~= 4.7px handles a
# single line of small tick-label text, verified empirically against a rendered report).
ROW_GRID_FRACTION = 0.86


def _row_geometry(row_idx: int, n: int) -> Tuple[str, str]:
    """`(top, height)` as ECharts percentage strings for grid row `row_idx` of `n`."""
    total = ROW_PX * n + EXTRA_PX + BOTTOM_PX
    top_pct = EXTRA_PX / total * 100
    bottom_pct = BOTTOM_PX / total * 100
    row_slot_pct = (100 - top_pct - bottom_pct) / n
    top = top_pct + row_idx * row_slot_pct
    height = row_slot_pct * ROW_GRID_FRACTION
    return f"{top:.4f}%", f"{height:.4f}%"


def _rgba(rgb: str, alpha: float) -> str:
    """`"rgb(r,g,b)"` -> `"rgba(r,g,b,alpha)"`."""
    return rgb.replace("rgb(", "rgba(").replace(")", f",{alpha})")


# --- KDE: the cross-language contract -------------------------------------------------
#
# GOLDEN VALUES (mirrored verbatim in `tests/test_plots_echarts.py::test_kde_golden_values`
# and, in Task 2.2, in a comment block at the top of
# `templates/echarts/src/js/plots/violin.js`). Any change to this formula must update both.
#
#   values = [1.0, 2.0, 3.0, 4.0, 5.0]
#   xs     = [1.0, 3.0, 5.0]
#   kde(values, xs) == [0.15916497933387785, 0.1802710624663249, 0.15916497933387785]
def kde(values: Sequence[float], xs: Sequence[float]) -> List[float]:
    """
    Epanechnikov-kernel density estimate with Silverman's rule-of-thumb bandwidth.

    Pure typed Python, no numpy: ported verbatim to JS in Task 2.2 (JS `kde()` must
    reproduce the golden values above bit-for-bit-close), and verbatim from
    `scripts/06_violin_final.py` lines 59-73. Keep this formula exactly as written; do
    not "clean it up" without re-deriving the golden values in both languages.
    """
    n = len(values)
    if n == 0:
        raise ValueError("kde: values must not be empty")
    mean = sum(values) / n
    variance = sum((v - mean) ** 2 for v in values) / n
    sd = variance**0.5
    bandwidth = max(1.06 * sd * n ** (-1 / 5), 1e-9)
    out = []
    for x in xs:
        density = 0.0
        for v in values:
            u = (x - v) / bandwidth
            if abs(u) <= 1:
                density += 0.75 * (1 - u * u)
        out.append(density / (n * bandwidth))
    return out


def _visible_metrics(dataset: Dataset) -> List[ColumnAnchor]:
    """
    Metrics to actually draw a row for: not hidden (a hidden table column stays hidden in
    the violin too, mirroring `Dataset.create_figure`'s `metrics = [m for m in self.metrics
    if not self.header_by_metric[m].hidden]`), and with at least one numeric value (a KDE
    needs numeric input; a metric that is entirely non-numeric strings has no violin row).
    """
    result = []
    for metric in dataset.metrics:
        header = dataset.header_by_metric[metric]
        if header.hidden:
            continue
        values = dataset.violin_value_by_sample_by_metric.get(metric, {})
        if any(isinstance(v, (int, float)) for v in values.values()):
            result.append(metric)
    return result


def _numeric_values(dataset: Dataset, metric: ColumnAnchor) -> List[float]:
    return [float(v) for v in dataset.violin_value_by_sample_by_metric[metric].values() if isinstance(v, (int, float))]


def _metric_range(header: ViolinColumn, values: List[float]) -> Tuple[float, float]:
    """
    Per-metric row's real x-axis range. Prefers the column's own configured axis range
    (`header.xaxis.range`, the exact `[dmin, dmax]` Plotly's per-row x-axis uses, see
    `Dataset.create` in `multiqc/plots/violin.py`), so a row's x-axis matches Plotly's
    exactly instead of always autoranging to the observed data spread (POLISH.md #10c).
    Falls back to a padded observed-value range (BUILD_PLAN.md Task 2.1), matching
    Plotly's own autorange, when the column has no configured range.
    """
    if header.xaxis.range is not None:
        lo, hi = float(header.xaxis.range[0]), float(header.xaxis.range[1])
        if hi > lo:
            return lo, hi

    lo, hi = min(values), max(values)
    span = hi - lo
    lo -= span * RANGE_PAD
    hi += span * RANGE_PAD
    if hi <= lo:
        # All values identical (span == 0): center a small synthetic window on the value
        # instead of pinning it to the row's left edge, so downstream (x - lo) / (hi - lo)
        # fractions (box/median placement) never divide by zero.
        lo, hi = values[0] - 0.5, values[0] + 0.5
    return lo, hi


def _metric_title(header: ViolinColumn) -> str:
    """
    Row label text for one metric (POLISH.md #7): titles/namespaces are authored with
    HTML entities (e.g. "Mapped &amp; paired") for Plotly's HTML-interpreting axis
    labels, but ECharts renders axis labels as plain text. `_clean_title_segment` (same
    helper `converter.py` uses for chart titles) unescapes entities and strips any stray
    tags.
    """
    if header.namespace:
        title = f"{header.namespace}: {header.title}"
    else:
        title = header.title
    return _clean_title_segment(title)


def _violin_polygon(values: List[float], lo: float, hi: float) -> List[List[float]]:
    """
    The KDE polygon for one row, in REAL x-space (each row draws on its own real-valued
    x-axis, POLISH.md #6) and a small y-space symmetric around 0 (the violin's density
    thickness; no row offset any more, since each row is its own grid): a closed shape
    formed by the density curve above the baseline followed by the same curve, mirrored,
    below it. Ported (minus the old 0..1 x-normalization) from `scripts/06_violin_final.py`
    lines 88-98.
    """
    span = hi - lo
    if max(values) == min(values):
        # Degenerate (zero-variance) metric: kde()'s Silverman bandwidth collapses to its
        # 1e-9 floor here, which would otherwise blow up into a single needle-thin,
        # full-height spike (POLISH.md #10d). Draw a flat row instead, matching Plotly's
        # own near-invisible rendering for constant-value metrics; the median tick drawn
        # over it (see series()) still marks the value.
        x = values[0]
        return [[x, 0.0], [x, 0.0]]

    xs = [lo + span * i / (N_BINS - 1) for i in range(N_BINS)]
    ys = kde(values, xs)
    ymax = max(ys) or 1.0
    top = [[x, (y / ymax) * ROW_HEIGHT] for x, y in zip(xs, ys)]
    bottom = [[x, -(y / ymax) * ROW_HEIGHT] for x, y in zip(reversed(xs), reversed(ys))]
    return top + bottom


def violin_polygons(dataset: Dataset) -> Dict[str, Dict[str, Any]]:
    """
    `{metric: {"poly": [[x, y], ...], "range": [lo, hi]}}`, computed from ALL samples in
    the dataset (not toolbox-filtered): the VIOLIN EXCEPTION described in the module
    docstring. `range` is the `[lo, hi]` real-value domain used for that metric's x-axis
    (the same pair `_violin_polygon` used to build `poly`'s real x coordinates).
    """
    result: Dict[str, Dict[str, Any]] = {}
    for metric in _visible_metrics(dataset):
        header = dataset.header_by_metric[metric]
        values = _numeric_values(dataset, metric)
        lo, hi = _metric_range(header, values)
        poly = _violin_polygon(values, lo, hi)
        result[str(metric)] = {"poly": poly, "range": [lo, hi]}
    return result


def _row_quartiles(values: List[float]) -> Tuple[float, float, float]:
    """
    Q1/median/Q3 for one metric row's raw values, via the same linear-interpolation
    quantile method `multiqc/plots/echarts/box.py` uses for box plots (reused directly,
    not duplicated), for the inner box + median line drawn over each violin
    (POLISH.md #10b).
    """
    sorted_values = sorted(values)
    return _quantile(sorted_values, 0.25), _quantile(sorted_values, 0.5), _quantile(sorted_values, 0.75)


def _violin_colors(header: ViolinColumn) -> Tuple[str, str, str]:
    """`[fill, stroke, box_fill]` for one metric row's violin: `fill` is the violin
    body's own semi-transparent color (`FILL_ALPHA`), `box_fill` is the inner Q1-Q3
    box's fill, the same stroke hue at a more solid alpha (`BOX_FILL_ALPHA`) so the box
    reads clearly against the lighter violin body."""
    if not header.color:
        return _DEFAULT_FILL, _DEFAULT_STROKE, _DEFAULT_BOX_FILL
    color = header.color
    # General-stats metric colors are bare "r,g,b" triples; color_to_rgb_string only
    # accepts rgb()/hex/named and would fall back to black, so wrap a bare triple first
    # (this matches normalizeColorToRGB on the JS side, keeping SSR colors in parity).
    if re.fullmatch(r"\s*\d+\s*,\s*\d+\s*,\s*\d+\s*", color):
        color = f"rgb({color})"
    stroke = color_to_rgb_string(color)  # "rgb(r,g,b)"
    fill = _rgba(stroke, FILL_ALPHA)
    box_fill = _rgba(stroke, BOX_FILL_ALPHA)
    return fill, stroke, box_fill


def _violin_render_series(
    poly: List[List[float]],
    q1: float,
    median: float,
    q3: float,
    row_idx: int,
    fill: str,
    stroke: str,
    box_fill: str,
) -> Dict[str, Any]:
    """
    One `custom` series per row whose `renderItem` draws the KDE polygon computed in
    Python, plus an inner Q1-Q3 box and median line (POLISH.md #10b), matching Plotly's
    `box_visible`/`meanline` violin config. `q1`/`median`/`q3` are real values (same
    x-space as `poly`); `xAxisIndex`/`yAxisIndex` bind this series to row `row_idx`'s own
    grid (POLISH.md #6/#8). This is the ONE place a real function needs to reach the SSR
    bundle: the `__FN__` sentinel body below is plain MultiQC-generated JS source (a
    JSON-encoded point list plus fixed numeric/color literals), never user data.
    """
    body = (
        f"var poly = {json.dumps(poly)};"
        "var pts = poly.map(function(p) { return api.coord(p); });"
        f"var bTL = api.coord([{json.dumps(q1)}, {json.dumps(-BOX_HALF_HEIGHT)}]);"
        f"var bBR = api.coord([{json.dumps(q3)}, {json.dumps(BOX_HALF_HEIGHT)}]);"
        f"var mTop = api.coord([{json.dumps(median)}, {json.dumps(-MEDIAN_HALF_HEIGHT)}]);"
        f"var mBot = api.coord([{json.dumps(median)}, {json.dumps(MEDIAN_HALF_HEIGHT)}]);"
        "return { type: 'group', children: ["
        "{ type: 'polygon', shape: { points: pts }, "
        f"style: {{ fill: {json.dumps(fill)}, stroke: {json.dumps(stroke)}, lineWidth: 1 }} }},"
        "{ type: 'rect', shape: { x: Math.min(bTL[0], bBR[0]), y: Math.min(bTL[1], bBR[1]), "
        "width: Math.abs(bBR[0] - bTL[0]), height: Math.abs(bBR[1] - bTL[1]) }, "
        f"style: {{ fill: {json.dumps(box_fill)}, stroke: {json.dumps(stroke)}, lineWidth: 1 }} }},"
        "{ type: 'line', shape: { x1: mTop[0], y1: mTop[1], x2: mBot[0], y2: mBot[1] }, "
        f"style: {{ stroke: {json.dumps(stroke)}, lineWidth: 2.5 }} }}"
        "] };"
    )
    return {
        "type": "custom",
        "coordinateSystem": "cartesian2d",
        "xAxisIndex": row_idx,
        "yAxisIndex": row_idx,
        "data": [0],
        "renderItem": {"__FN__": True, "body": body},
        "silent": True,
        "z": 1,
    }


def _annotation_render_series(row_idx: int, lo_obs: float, hi_obs: float) -> Dict[str, Any]:
    """
    One `custom` series per row drawing the observed min/max text labels at the violin's
    actual real-valued extent (`lo_obs`/`hi_obs`), vertically at the row's center (`y=0`).
    Plotly still shows these value labels at the violin ends even though its own axis
    already shows real values, so they're kept here too (POLISH.md #6 note).
    """
    lo_label = json.dumps(f"{lo_obs:g}")
    hi_label = json.dumps(f"{hi_obs:g}")
    body = (
        f"var L = api.coord([{json.dumps(lo_obs)}, 0]); var R = api.coord([{json.dumps(hi_obs)}, 0]);"
        "return { type: 'group', children: ["
        "{ type: 'text', style: { text: " + lo_label + ", x: L[0] - 6, y: L[1], "
        "textAlign: 'right', textVerticalAlign: 'middle', fontSize: 10, fill: '#888' } }, "
        "{ type: 'text', style: { text: " + hi_label + ", x: R[0] + 6, y: R[1], "
        "textAlign: 'left', textVerticalAlign: 'middle', fontSize: 10, fill: '#888' } }"
        "] };"
    )
    return {
        "type": "custom",
        "coordinateSystem": "cartesian2d",
        "xAxisIndex": row_idx,
        "yAxisIndex": row_idx,
        "data": [0],
        "renderItem": {"__FN__": True, "body": body},
        "silent": True,
        "z": 3,
    }


def _scatter_series_for_metric(dataset: Dataset, metric: ColumnAnchor, row_idx: int) -> Dict[str, Any]:
    """
    ONE beeswarm `scatter` series for this metric (not one shared series across all
    metrics, unlike the reference script): `jitter`/`jitterOverlap` gives ECharts 6's
    native beeswarm spread. Items are `{"value": [real_value, 0], "name": sample}` (real
    value, not normalized; POLISH.md #6); per-sample coloring and the interactive
    tooltip are added by the interactive JS builder in Task 2.2, not here.
    """
    header = dataset.header_by_metric[metric]
    data: List[Dict[str, Any]] = []
    if header.show_points:
        if header.show_only_outliers:
            values_by_sample = dataset.scatter_value_by_sample_by_metric.get(metric, {})
        else:
            values_by_sample = dataset.violin_value_by_sample_by_metric.get(metric, {})
        for sample, value in values_by_sample.items():
            if not isinstance(value, (int, float)):
                continue
            data.append({"value": [float(value), 0], "name": str(sample)})

    return {
        "type": "scatter",
        "name": str(metric),
        "xAxisIndex": row_idx,
        "yAxisIndex": row_idx,
        "data": data,
        "symbolSize": 6,
        "jitter": 22,
        "jitterOverlap": False,
        "itemStyle": {"color": _SCATTER_COLOR, "borderColor": "#fff", "borderWidth": 0.5},
        "z": 2,
    }


# Rough average glyph width (px) for the row-title label font (fontSize 12, default
# sans-serif), used to reserve a left margin generous enough to avoid clipping.
_TITLE_CHAR_PX = 6.6
_MIN_GRID_LEFT = 60
_MAX_GRID_LEFT = 280


def _grid_left(titles: List[str]) -> int:
    """
    Shared left inset (px) for every row's grid, reserved for the row title label
    (`_row_axes`'s yAxis), sized from the LONGEST visible title so every row's plot area
    starts at the same x (a column of labels reads better aligned than staggered).
    `containLabel: true` can still grow this further for an unusually long title, but
    is not trusted as the ONLY source of this margin: `containLabel`'s automatic text
    measurement runs on real Canvas metrics in a browser (accurate) but on a crude
    fallback estimate in the SSR path's headless MiniRacer/V8 context (no real Canvas
    API), which can under-measure a long title enough to clip it off the left edge of
    the exported PNG; computing this floor from the title's own character count keeps
    the SSR and interactive paths both correct rather than only the browser one.
    """
    if not titles:
        return _MIN_GRID_LEFT
    max_len = max(len(t) for t in titles)
    return int(min(_MAX_GRID_LEFT, max(_MIN_GRID_LEFT, max_len * _TITLE_CHAR_PX + 14)))


def _row_axes(
    row_idx: int, n: int, header: ViolinColumn, lo: float, hi: float, title: str, grid_left: int
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """
    `(grid, xAxis, yAxis)` for one metric row: a real-valued, SI-formatted x-axis
    (POLISH.md #6) and a hidden value y-axis whose single tick (at `y=0`, the exact
    vertical center of the row) carries the metric title (POLISH.md #1, the VALUE-AXIS
    TRICK described in the module docstring, now one instance per row instead of one
    shared instance for every row). `grid_left` (see `_grid_left`) is shared across every
    row so the plot areas all start at the same x.
    """
    top, height = _row_geometry(row_idx, n)
    grid = {"top": top, "height": height, "left": grid_left, "right": 16, "containLabel": True}

    suffix = str(header.xaxis.ticksuffix) if header.xaxis.ticksuffix else ""
    x_axis = {
        "type": "value",
        "gridIndex": row_idx,
        "min": lo,
        "max": hi,
        "axisLabel": {"fontSize": 10, "formatter": {"__FN__": True, "body": _si_axis_formatter_body(suffix)}},
        "axisLine": {"show": True},
        "axisTick": {"show": True},
        "splitLine": {"show": False},
    }
    y_axis = {
        "type": "value",
        "gridIndex": row_idx,
        "min": -0.5,
        "max": 0.5,
        "interval": 0.5,
        "axisLabel": {
            "fontSize": 12,
            "align": "right",
            "verticalAlign": "middle",
            "formatter": {"__FN__": True, "body": f"return Math.abs(v) < 1e-6 ? {json.dumps(title)} : '';"},
        },
        "axisTick": {"show": False},
        "axisLine": {"show": False},
        "splitLine": {"show": False},
    }
    return grid, x_axis, y_axis


def layout_option(plot: "Plot[Any, Any]", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one violin dataset, minus the KDE/annotation/scatter
    `series` (built by `series()` below for the SSR/get_option path). `grid`/`xAxis`/
    `yAxis` are each a LIST with one entry per visible metric (PLOTLY-STYLE PER-ROW
    SUBPLOTS), not the single shared grid/axis pair every other type gets from
    `convert_layout`. See the module docstring for why the yAxis entries carry a `__FN__`
    formatter sentinel (the VIOLIN EXCEPTION).
    """
    option = convert_layout(plot.layout, dataset.layout)

    metrics = _visible_metrics(dataset)
    n = len(metrics)
    titles = [_metric_title(dataset.header_by_metric[m]) for m in metrics]
    grid_left = _grid_left(titles)

    grids: List[Dict[str, Any]] = []
    x_axes: List[Dict[str, Any]] = []
    y_axes: List[Dict[str, Any]] = []
    for row_idx, metric in enumerate(metrics):
        header = dataset.header_by_metric[metric]
        values = _numeric_values(dataset, metric)
        lo, hi = _metric_range(header, values)
        grid, x_axis, y_axis = _row_axes(row_idx, n, header, lo, hi, titles[row_idx], grid_left)
        grids.append(grid)
        x_axes.append(x_axis)
        y_axes.append(y_axis)

    option["grid"] = grids
    option["xAxis"] = x_axes
    option["yAxis"] = y_axes
    option["tooltip"]["trigger"] = "item"

    # No cross spike lines (POLISH.md #13): the y-axis in every row is still the
    # VALUE-AXIS TRICK (module docstring), a fake symmetric thickness scale, not real
    # data (only x is real now); a value label following the cursor on that axis would
    # show numbers that mean nothing to the viewer. The crosshair mouse *cursor* (set
    # globally in JS, see echarts-plotting.js's buildCurrentOption) still applies to
    # violin like every other type; only the axis-tracking guide lines are skipped here.
    option["tooltip"]["axisPointer"] = {"show": False}
    option["axisPointer"] = {"show": False}

    # Plotly-style click+drag box-zoom, one real x-axis per row (POLISH.md #8): the
    # toolbox dataZoom feature spans every row's xAxisIndex, so a drag-select inside any
    # row's grid zooms only THAT row's x-axis; one `inside` dataZoom per row holds each
    # row's own current zoom range (mirrors `converter.zoom_option`'s "one inside
    # dataZoom per zoomable axis" reasoning, generalized from 1 axis to N). No yAxisIndex
    # anywhere: the y-axis is the fake per-row thickness scale, never worth zooming.
    option["toolbox"] = {
        "show": True,
        "top": "150%",
        "feature": {"dataZoom": {"show": True, "xAxisIndex": list(range(n)), "yAxisIndex": []}},
    }
    option["dataZoom"] = [
        {"type": "inside", "xAxisIndex": [i], "zoomOnMouseWheel": False, "moveOnMouseWheel": False} for i in range(n)
    ]

    return option


def series(dataset: Dataset, pconfig: TableConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    `n_metric` violin (KDE polygon + inner box + median) `custom` series, then `n_metric`
    min/max-annotation `custom` series, then `n_metric` beeswarm `scatter` series (one per
    metric; see `_scatter_series_for_metric`), each bound to its own row's
    `xAxisIndex`/`yAxisIndex`. This is the SSR/get_option (non-toolbox) path; the
    interactive path is `EchartsViolinPlot.buildSeries()`
    (`templates/echarts/src/js/plots/violin.js`, Task 2.2).

    `pconfig`/`is_pct` are accepted for dispatch-signature parity with `bar.series`; violin
    plots have no percentage switch, so `is_pct` is unused here.
    """
    metrics = _visible_metrics(dataset)

    violin_series = []
    annotation_series = []
    scatter_series = []
    for row_idx, metric in enumerate(metrics):
        header = dataset.header_by_metric[metric]
        values = _numeric_values(dataset, metric)
        lo, hi = _metric_range(header, values)

        poly = _violin_polygon(values, lo, hi)
        fill, stroke, box_fill = _violin_colors(header)
        q1, median, q3 = _row_quartiles(values)
        violin_series.append(_violin_render_series(poly, q1, median, q3, row_idx, fill, stroke, box_fill))

        lo_obs, hi_obs = min(values), max(values)
        annotation_series.append(_annotation_render_series(row_idx, lo_obs, hi_obs))

        scatter_series.append(_scatter_series_for_metric(dataset, metric, row_idx))

    return violin_series + annotation_series + scatter_series


def axis_data(dataset: Dataset, pconfig: TableConfig) -> Optional[List[Tuple[str, List[str]]]]:
    """
    Always `None`: every row's x/y axis is `type: "value"` (a real-valued x-axis, a
    VALUE-AXIS TRICK y-axis), so there is no category `data` array for the toolbox to
    fill in.
    """
    return None


def mark_count(dataset: Dataset) -> int:
    """Total scatter (beeswarm) points across all visible metrics."""
    total = 0
    for metric in _visible_metrics(dataset):
        header = dataset.header_by_metric[metric]
        if not header.show_points:
            continue
        values_by_sample = (
            dataset.scatter_value_by_sample_by_metric.get(metric, {})
            if header.show_only_outliers
            else dataset.violin_value_by_sample_by_metric.get(metric, {})
        )
        total += sum(1 for v in values_by_sample.values() if isinstance(v, (int, float)))
    return total
