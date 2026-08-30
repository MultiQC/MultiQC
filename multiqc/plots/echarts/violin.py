"""
ECharts option builder for `PlotType.VIOLIN` (`multiqc/plots/violin.py`).

Ports `multiqc-echarts-exploration/scripts/06_violin_final.py` faithfully, then extends
it with PLOTLY-STYLE PER-ROW SUBPLOTS: one ECharts `grid` (with its own `xAxis`/`yAxis`
pair) per visible metric, stacked vertically, matching `multiqc/plots/violin.py`'s own
`Dataset.create_figure` (one Plotly subplot per metric, real x-axis values, compact
per-row height). Each row is drawn as a hand-built KDE polygon via a `custom` series
`renderItem`, plus a beeswarm `scatter` series with v6's native `jitter`. (An earlier
version also drew a third `custom` series per row for min/max text labels next to the
violin ends; removed, since each row already has its own real-valued x-axis showing that
same range, so the labels were pure duplication and read as if the beeswarm dots
themselves carried on-canvas text.)

Per-row geometry: `_row_geometry` mirrors the same formula duplicated in
`templates/echarts/src/js/plots/violin.js` (`rowGeometry`), so a row's grid percentage
split is identical whether computed here (SSR/`get_option`) or in the browser
(interactive). Because grid percentages are resolved against whatever pixel height the
chart is actually given, this self-scales correctly for both the SSR/flat path (container
height comes from the shared Plotly-style `plot.layout.height`, unrelated to this file)
and the interactive path (the JS renderer sets a compact, row-count-scaled CSS height
itself; see `EchartsViolinPlot.buildSeries()`).

Since each row now gets its OWN real-valued x-axis, a row's KDE polygon/box/median/
scatter coordinates are drawn directly in REAL x (no more per-row 0..1
normalization); only the y coordinate stays a small symmetric offset around 0 (the
violin's density thickness), identical in scale to the old `ROW_HEIGHT`/`BOX_HALF_HEIGHT`/
`MEDIAN_HALF_HEIGHT` constants, just no longer shifted by a row index.

Two deliberate exceptions to the general "ECharts model->JSON contract"
(`multiqc-echarts-exploration/BUILD_PLAN.md`), both called out there by name:

1. ROW TITLE AS A FIXED `title` COMPONENT (formerly a "VALUE-AXIS TRICK" that drew the
   title as the yAxis tick label at y=0): that approach made `containLabel: true` grow
   each row's `grid.left` by however much THAT row's own title text measured, so rows
   started at a different x with a different width (measured 329..454px x, 486..599px
   width across a real fastqc+samtools report). The yAxis is now a plain hidden `type:
   "value"` axis spanning a fixed `[-0.5, 0.5]` (`axisLabel.show: false`, contributing
   NOTHING to the grid's inset), purely so `renderItem`s below have a small vertical
   coordinate range to draw the violin's thickness against. The title itself is drawn by
   `_title_option`, ONE MORE entry in ECharts' native `title` component array (`option
   ["title"]`, alongside the chart's own main title from `convert_layout`): `left`/`top`
   position it at a FIXED pixel/percentage independent of any grid or axis, with
   `textAlign: "right"`/`textVerticalAlign: "middle"` right-aligning and vertically
   centering it (POLISH.md #1 still holds). This was tried first as a `custom` series
   `renderItem` (reading `params.coordSys`, the grid's own pixel rect, to stay fixed
   regardless of zoom): that rendered correctly in a live browser, but ECharts' SSR mode
   (`static_export.py`, `ssr: true`) silently drops ANY manually-drawn `type: "text"`
   element, whether from a `custom` series `renderItem` OR the `graphic` component,
   because SSR has no real Canvas 2D context to measure arbitrary text with; only text
   from ECharts' own built-in components (axis labels, `title`, `legend`, ...) survives
   SSR, since those ship their own fallback measurement. `title` array entries are
   exactly such a built-in component, so this is the one placement that renders
   identically in both the browser and the MiniRacer/resvg SSR path (verified by
   rendering the real fastqc+samtools general-stats violin through both). Every grid
   still shares the same `left`/`right` and `containLabel: false` (see
   `_grid_left`/`_row_axes`), so the pixel position `_title_option` computes for each row
   lines up with that row's grid exactly, in a single column instead of staggering.

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

# Fill opacity for the violin body (POLISH.md #10a). Measured directly from Plotly's own
# rendered general-stats violin (`gsDiv.data[i].fillcolor`, live browser inspection, not
# guessed): Plotly's light-mode fillcolor is `rgba(<header.color or #999999>, 0.5)` (see
# `templates/default/src/js/plots/violin.js`'s client-side `isDarkMode ? 0.3 : 0.5`; the
# report always renders in light mode there, so 0.5 is the value actually shown, not 0.3
# as an earlier version of this comment assumed). The SSR path here always renders
# against the light theme too (see module docstring), so both this module and the
# interactive JS path fix on the same 0.5 value.
FILL_ALPHA = 0.5

# Default violin body color (POLISH.md #5): matches Plotly's own general-stats violin
# default, `fillcolor="#999999"` (see `multiqc/plots/violin.py`'s `Dataset.create`), kept
# as rgb(...) so _rgba() below applies uniformly. Only used when a metric has no
# configured `header.color`; a real per-metric color always takes priority (see
# `_violin_colors` below).
_DEFAULT_STROKE = "rgb(153,153,153)"
_DEFAULT_FILL = f"rgba(153,153,153,{FILL_ALPHA})"

# Beeswarm point defaults (POLISH.md #9/FIX #9): matches `multiqc/plots/violin.py`'s own
# `Dataset.create` scatter_trace_params exactly (`marker: {size: 4, color: ..., opacity: 1}`,
# no marker line/border), so the echarts beeswarm reads the same as Plotly's. `_SCATTER_SIZE`
# is the non-highlighted point diameter; the SSR path never renders a highlighted state
# (that's an interactive-only toolbox feature), so unlike the JS port there is only one size.
_SCATTER_SIZE = 4
# Plotly's own choice: black when any metric has a custom color (so grey dots don't clash
# with a color-coded violin), else blue (so a plain grey violin's dots read as interactive).
# The SSR/flat path always renders against the registered "multiqc-light" theme (see
# static_export.py), so unlike the interactive JS port (which also handles a dark viewer
# theme) this never needs the dark-mode swap Plotly's own JS applies client-side.
_SCATTER_COLOR_COLORED = "#000000"
_SCATTER_COLOR_PLAIN = "#0b79e6"

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
# Fraction of a row's slot handed to ECharts as the grid box (`containLabel: false`, see
# `_row_axes`: the grid rect is exactly this box, not auto-shrunk, so the row's own
# x-axis tick labels are drawn just BELOW it, in the remaining 1 - 0.86 of the slot, as a
# small gap before the next row's grid starts; 34px * (1 - 0.86) ~= 4.7px handles a
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


def _row_center_pct(row_idx: int, n: int) -> float:
    """
    Vertical center of grid row `row_idx` of `n`, as a float percentage of the container
    height (same basis as `_row_geometry`'s `top`/`height`, duplicated here rather than
    parsed back out of those formatted strings, matching this module's usual
    cross-language-duplication style, e.g. `kde()`): the `title` component drawing that
    row's metric label (`_title_option`) is positioned here so it stays vertically
    centered on the row regardless of row count (POLISH.md #1).
    """
    total = ROW_PX * n + EXTRA_PX + BOTTOM_PX
    top_pct = EXTRA_PX / total * 100
    bottom_pct = BOTTOM_PX / total * 100
    row_slot_pct = (100 - top_pct - bottom_pct) / n
    row_top_pct = top_pct + row_idx * row_slot_pct
    row_height_pct = row_slot_pct * ROW_GRID_FRACTION
    return row_top_pct + row_height_pct / 2


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


def _violin_colors(header: ViolinColumn) -> Tuple[str, str]:
    """`[fill, stroke]` for one metric row's violin: `fill` is the violin body's own
    semi-transparent color (`FILL_ALPHA`), used for the body AND the inner Q1-Q3 box
    (Plotly's own violin draws its box with the exact same `fillcolor` as the body, see
    the module docstring); `stroke` (the solid hue, no alpha) is used only for the
    median line, never as a border (Plotly's violin outline/box/meanline all render with
    `line.width: 0`, i.e. no border at all, measured directly from its rendered
    `_fullData`; a bordered polygon reads far more saturated/"brighter" than Plotly's
    soft alpha-only fill, see `_violin_render_series`)."""
    if not header.color:
        return _DEFAULT_FILL, _DEFAULT_STROKE
    color = header.color
    # General-stats metric colors are bare "r,g,b" triples; color_to_rgb_string only
    # accepts rgb()/hex/named and would fall back to black, so wrap a bare triple first
    # (this matches normalizeColorToRGB on the JS side, keeping SSR colors in parity).
    if re.fullmatch(r"\s*\d+\s*,\s*\d+\s*,\s*\d+\s*", color):
        color = f"rgb({color})"
    stroke = color_to_rgb_string(color)  # "rgb(r,g,b)"
    fill = _rgba(stroke, FILL_ALPHA)
    return fill, stroke


def _violin_render_series(
    poly: List[List[float]],
    q1: float,
    median: float,
    q3: float,
    row_idx: int,
    fill: str,
    stroke: str,
) -> Dict[str, Any]:
    """
    One `custom` series per row whose `renderItem` draws the KDE polygon computed in
    Python, plus an inner Q1-Q3 box and median line (POLISH.md #10b), matching Plotly's
    `box_visible`/`meanline` violin config. `q1`/`median`/`q3` are real values (same
    x-space as `poly`); `xAxisIndex`/`yAxisIndex` bind this series to row `row_idx`'s own
    grid (POLISH.md #6/#8). This is the ONE place a real function needs to reach the SSR
    bundle: the `__FN__` sentinel body below is plain MultiQC-generated JS source (a
    JSON-encoded point list plus fixed numeric/color literals), never user data.

    Neither the polygon nor the box has a `stroke` (matches Plotly's measured
    `line.width: 0`, see `_violin_colors`): the box reuses the SAME `fill` as the body
    (also matching Plotly, whose box `fillcolor` is identical to the trace's own), so it
    reads as a slightly denser patch purely from the two semi-transparent fills
    overlapping, not from an added border. Only the median line keeps a solid `stroke`
    (Plotly's own `meanline.width` measures 0, i.e. invisible; MultiQC keeps a thin
    visible tick here since a fully invisible median would be a regression, not a match).
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
        f"style: {{ fill: {json.dumps(fill)} }} }},"
        "{ type: 'rect', shape: { x: Math.min(bTL[0], bBR[0]), y: Math.min(bTL[1], bBR[1]), "
        "width: Math.abs(bBR[0] - bTL[0]), height: Math.abs(bBR[1] - bTL[1]) }, "
        f"style: {{ fill: {json.dumps(fill)} }} }},"
        "{ type: 'line', shape: { x1: mTop[0], y1: mTop[1], x2: mBot[0], y2: mBot[1] }, "
        f"style: {{ stroke: {json.dumps(stroke)}, lineWidth: 2 }} }}"
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


def _scatter_base_color(dataset: Dataset) -> str:
    """
    One base beeswarm color for the whole dataset (mirrors `multiqc/plots/violin.py`'s
    `Dataset.create`: `"#000000" if any(h.color is not None for h in
    ds.header_by_metric.values()) else "#0b79e6"`), so a color-coded violin's dots don't
    visually clash with it, while a plain grey violin gets a blue "this is interactive"
    hint. Computed once per dataset, not per row: every row's beeswarm shares it.
    """
    if any(h.color is not None for h in dataset.header_by_metric.values()):
        return _SCATTER_COLOR_COLORED
    return _SCATTER_COLOR_PLAIN


def _scatter_series_for_metric(dataset: Dataset, metric: ColumnAnchor, row_idx: int, base_color: str) -> Dict[str, Any]:
    """
    ONE beeswarm `scatter` series for this metric (not one shared series across all
    metrics, unlike the reference script): `jitter`/`jitterOverlap` gives ECharts 6's
    native beeswarm spread. Items are `{"value": [real_value, 0], "name": sample}` (real
    value, not normalized; POLISH.md #6); per-sample coloring and the interactive
    tooltip are added by the interactive JS builder in Task 2.2, not here.

    Point style (FIX #9) matches Plotly's own beeswarm marker exactly: size 4, solid
    (`opacity: 1`), no border (Plotly's `scatter_trace_params` never sets a marker line).
    `label.show: False` (FIX #2) is explicit, not just relying on ECharts' own default,
    so no per-point text ever draws next to a dot outside the hover tooltip. `z: 3` is
    ABOVE `_violin_render_series`'s `z: 1` (the KDE polygon/box/median), so points always
    draw on top of the (now un-bordered, semi-transparent) violin body, never obscured by
    it.
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
        "symbolSize": _SCATTER_SIZE,
        "jitter": 22,
        "jitterOverlap": False,
        "itemStyle": {"color": base_color, "opacity": 1},
        "label": {"show": False},
        "z": 3,
    }


# Rough average glyph width (px) for the row-title label font (fontSize 12, default
# sans-serif), used to reserve a left margin generous enough to avoid clipping. Measured
# via canvas.measureText() against real MultiQC general-stats titles (~5.3px/char for the
# longest title in a fastqc+samtools report); 6.0 keeps a safety margin above that
# without reserving the excessive empty space the old 6.6 did (ROW ALIGNMENT fix).
_TITLE_CHAR_PX = 6.0
_MIN_GRID_LEFT = 60
_MAX_GRID_LEFT = 280
# Gap (px) between the right-aligned title text and the row's plot area (ROW ALIGNMENT
# fix): the title is drawn `_TITLE_GUTTER_PAD` px left of the grid's own rect (see
# `_title_option`). Wide enough to clear the row's OWN x-axis first tick label
# (e.g. "0%"/"0 M"/"627 bp"), which straddles that same left edge (a value axis has no
# `boundaryGap`) and would otherwise overlap the title text now that `containLabel`/
# `outerBoundsMode` no longer reserve room for it (measured worst case ~30px wide at
# fontSize 10, so half its width plus a little slack comfortably fits in 24px).
_TITLE_GUTTER_PAD = 24
# Small slack (px) between the title text's own estimated width and where
# `_TITLE_GUTTER_PAD` starts, absorbing `_TITLE_CHAR_PX`'s per-glyph estimation error.
_TITLE_TEXT_SLACK = 4
# Mirrors `theme.py`'s `_TEXT_COLOR` (the SSR theme's axisLabel/title color): the SSR
# path always renders against the light theme (see that module's own docstring), so this
# is a plain literal rather than a cross-module import of a private theme.py constant.
_TITLE_TEXT_COLOR = "#333333"

# Per-row x-axis tick label color, measured directly from Plotly's own rendered general
# stats violin (`gsDiv.layout.xaxis.tickfont.color`, live browser inspection): unlike
# `_TITLE_TEXT_COLOR` (the row's own name, deliberately dark/prominent), Plotly's axis
# tick text is a subdued mid-grey. The SSR theme's own `valueAxis.axisLabel.color`
# (`theme.py`'s `_TEXT_COLOR`, "#333333") is darker/more prominent than this, which is
# part of why the SSR violin's axis previously read brighter/heavier than Plotly's.
_AXIS_LABEL_COLOR = "rgba(80,80,80,1)"


def _grid_left(titles: List[str]) -> int:
    """
    Shared left inset (px) for every row's grid, reserved for the row title text (drawn
    by `_title_option`, not the yAxis any more, see the module docstring's point 1),
    sized from the LONGEST visible title so every row's plot area starts at the same x (a
    column of labels reads better aligned than staggered). Computed from the title's own
    character count (not measured) so the SSR path's headless MiniRacer/V8 context (no
    real Canvas API for accurate text metrics) is just as correct as the browser.
    """
    if not titles:
        return _MIN_GRID_LEFT
    max_len = max(len(t) for t in titles)
    reserved = max_len * _TITLE_CHAR_PX + _TITLE_TEXT_SLACK + _TITLE_GUTTER_PAD
    return int(min(_MAX_GRID_LEFT, max(_MIN_GRID_LEFT, reserved)))


def _row_axes(
    row_idx: int, n: int, header: ViolinColumn, lo: float, hi: float, grid_left: int
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """
    `(grid, xAxis, yAxis)` for one metric row: a real-valued, SI-formatted x-axis
    (POLISH.md #6) and a fully hidden value y-axis (just a `[-0.5, 0.5]` coordinate range
    for `renderItem`s to draw the violin's thickness against; the title text is drawn
    separately, see `_title_option`). `grid_left` (see `_grid_left`) and `right`
    are IDENTICAL for every row, and `containLabel: False` means nothing (no more visible
    yAxis label, and x-axis tick labels are short/fixed-position) can grow one row's rect
    differently from another's: every row's `coordinateSystem.getRect()` is therefore
    equal in both `x` and `width` (the ROW ALIGNMENT fix; see the module docstring).
    """
    top, height = _row_geometry(row_idx, n)
    grid = {
        "top": top,
        "height": height,
        "left": grid_left,
        "right": 16,
        "containLabel": False,
        # ECharts 6 also auto-nudges a grid's rect inward, independently of
        # `containLabel`, to keep axis labels from overflowing the chart's OUTER
        # bounds (`grid.outerBoundsMode` default `"auto"`, see `GridModel.js`'s
        # `OUTER_BOUNDS_DEFAULT`): since a row's x-axis tick label LENGTH varies by
        # metric (e.g. "20%" vs "18.39 M"), this on its own reproduces the exact same
        # per-row rect inequality bug `containLabel` caused, just via a different
        # ECharts 6 feature. `"none"` disables it, verified empirically (see the ROW
        # ALIGNMENT fix) to make every row's rect byte-for-byte identical.
        "outerBoundsMode": "none",
    }

    suffix = str(header.xaxis.ticksuffix) if header.xaxis.ticksuffix else ""
    x_axis = {
        "type": "value",
        "gridIndex": row_idx,
        "min": lo,
        "max": hi,
        "axisLabel": {
            "fontSize": 10,
            "color": _AXIS_LABEL_COLOR,
            "formatter": {"__FN__": True, "body": _si_axis_formatter_body(suffix)},
        },
        # No axis line/ticks (measured Plotly's own violin x-axis: `showline: false`,
        # `ticks: ""`; only the tick label text is shown). A full-width axis line PLUS
        # tick marks under every row, at any color, reads far more present/"brighter"
        # than Plotly's bare label text, and was the main source of the axis-brightness
        # mismatch (see the module docstring/PR notes).
        "axisLine": {"show": False},
        "axisTick": {"show": False},
        "splitLine": {"show": False},
    }
    y_axis = {
        "type": "value",
        "gridIndex": row_idx,
        "min": -0.5,
        "max": 0.5,
        "interval": 0.5,
        "axisLabel": {"show": False},
        "axisTick": {"show": False},
        "axisLine": {"show": False},
        "splitLine": {"show": False},
    }
    return grid, x_axis, y_axis


def _title_option(row_idx: int, n: int, title: str, grid_left: int) -> Dict[str, Any]:
    """
    One entry in ECharts' native `title` component array (`option["title"]`) drawing the
    row's metric TITLE, fixed at `grid_left - _TITLE_GUTTER_PAD` px from the container's
    left edge and vertically centered on the row (`_row_center_pct`), right-aligned so its
    text ends exactly there. `padding: 0`/`fontWeight: "normal"` strip the `title`
    component's own defaults (a title normally reserves internal padding and renders
    bold), so `left`/`top` land pixel-exact and the row title reads like a plain label,
    not a heading.

    A NATIVE ECharts component, not a `custom` series `renderItem` or `graphic` element:
    tried both of those first, since either can read a fixed pixel position immune to
    zoom (see the module docstring's point 1), but ECharts' SSR mode
    (`static_export.py`) silently drops manually-drawn `type: "text"` elements from
    EITHER of those, because it has no real Canvas 2D context to measure arbitrary text
    with. `title` array entries are a built-in component (same family as the chart's own
    main title, axis labels, and legend) that ships its own SSR-safe measurement
    fallback, so this is the one placement verified to render in both the browser and the
    MiniRacer/resvg SSR path.
    """
    center_pct = _row_center_pct(row_idx, n)
    return {
        "text": title,
        "left": grid_left - _TITLE_GUTTER_PAD,
        "top": f"{center_pct:.4f}%",
        "textAlign": "right",
        "textVerticalAlign": "middle",
        "padding": 0,
        "textStyle": {"fontSize": 12, "fontWeight": "normal", "color": _TITLE_TEXT_COLOR},
    }


def layout_option(plot: "Plot[Any, Any]", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one violin dataset, minus the KDE/title/scatter
    `series` (built by `series()` below for the SSR/get_option path). `grid`/
    `xAxis`/`yAxis` are each a LIST with one entry per visible metric (PLOTLY-STYLE
    PER-ROW SUBPLOTS), not the single shared grid/axis pair every other type gets from
    `convert_layout`. See the module docstring for why the xAxis entries carry a `__FN__`
    formatter sentinel (the VIOLIN EXCEPTION); the yAxis is a plain hidden axis (no
    formatter, no label at all, see `_row_axes`).
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
        grid, x_axis, y_axis = _row_axes(row_idx, n, header, lo, hi, grid_left)
        grids.append(grid)
        x_axes.append(x_axis)
        y_axes.append(y_axis)

    option["grid"] = grids
    option["xAxis"] = x_axes
    option["yAxis"] = y_axes
    option["tooltip"]["trigger"] = "item"

    # Row titles (module docstring point 1): ECharts' `title` option accepts either one
    # component or a LIST of them, so the row titles ride alongside the chart's own main
    # title (set by `convert_layout` above) rather than replacing it.
    main_title = option.get("title")
    option["title"] = ([main_title] if main_title else []) + [
        _title_option(row_idx, n, titles[row_idx], grid_left) for row_idx in range(n)
    ]

    # No cross spike lines (POLISH.md #13): the y-axis in every row is a fake symmetric
    # thickness scale, not real data (only x is real, see the module docstring); a value
    # label following the cursor on that axis would
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
        {
            "type": "inside",
            "xAxisIndex": [i],
            "zoomOnMouseWheel": False,
            "moveOnMouseWheel": False,
            # `filterMode: "none"` (default is "filter"): ECharts' default dataZoom
            # behavior REMOVES a series' data points once their x-value falls outside the
            # zoomed range, which for `_violin_render_series` (a dummy `data: [0]`, not a
            # real x-coordinate; its renderItem draws from closure state, not from that
            # data value) makes ECharts drop the whole row's polygon the moment a zoom
            # narrows the axis so 0 falls outside it. "none" zooms the axis view without
            # filtering any series data, fixing a pre-existing bug where zooming any row
            # already made its violin polygon vanish this way.
            "filterMode": "none",
        }
        for i in range(n)
    ]

    return option


def series(dataset: Dataset, pconfig: TableConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    `n_metric` violin (KDE polygon + inner box + median) `custom` series, then `n_metric`
    beeswarm `scatter` series (one per metric; see `_scatter_series_for_metric`), each
    bound to its own row's `xAxisIndex`/`yAxisIndex`. This is the SSR/get_option
    (non-toolbox) path; the interactive path is `EchartsViolinPlot.buildSeries()`
    (`templates/echarts/src/js/plots/violin.js`, Task 2.2). The row TITLES are not a
    series at all (see `_title_option`/`layout_option`), so they are not built here.

    There is deliberately no min/max-annotation series any more: each row already has its
    own real-valued x-axis (POLISH.md #6), so a text label repeating the same min/max
    values right next to the dots was pure duplication, and read as if the beeswarm points
    themselves carried on-canvas labels (they never did; that was this now-removed series).

    `pconfig`/`is_pct` are accepted for dispatch-signature parity with `bar.series`; violin
    plots have no percentage switch, so `is_pct` is unused here.
    """
    metrics = _visible_metrics(dataset)
    base_color = _scatter_base_color(dataset)

    violin_series = []
    scatter_series = []
    for row_idx, metric in enumerate(metrics):
        header = dataset.header_by_metric[metric]
        values = _numeric_values(dataset, metric)
        lo, hi = _metric_range(header, values)

        poly = _violin_polygon(values, lo, hi)
        fill, stroke = _violin_colors(header)
        q1, median, q3 = _row_quartiles(values)
        violin_series.append(_violin_render_series(poly, q1, median, q3, row_idx, fill, stroke))

        scatter_series.append(_scatter_series_for_metric(dataset, metric, row_idx, base_color))

    return violin_series + scatter_series


def axis_data(dataset: Dataset, pconfig: TableConfig) -> Optional[List[Tuple[str, List[str]]]]:
    """
    Always `None`: every row's x/y axis is `type: "value"` (a real-valued x-axis, a
    fully hidden y-axis), so there is no category `data` array for the toolbox to fill in.
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
