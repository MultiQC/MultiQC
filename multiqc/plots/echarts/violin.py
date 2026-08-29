"""
ECharts option builder for `PlotType.VIOLIN` (`multiqc/plots/violin.py`).

Ports `multiqc-echarts-exploration/scripts/06_violin_final.py` faithfully. Violin is the
one plot type ECharts has no native chart for: each row is drawn as a hand-built polygon
via a `custom` series `renderItem`, plus a beeswarm `scatter` series with v6's native
`jitter`, plus a `custom` series per row for the min/max text annotations.

Two deliberate exceptions to the general "ECharts model->JSON contract"
(`multiqc-echarts-exploration/BUILD_PLAN.md`), both called out there by name:

1. VALUE-AXIS TRICK: `yAxis` is a `type: "value"` axis (not `"category"`), because a
   category axis snaps the violins' fractional row offsets (the KDE polygon extends
   above/below its row's integer index) to the nearest integer tick. `layout_option`
   therefore emits `yAxis.axisLabel.formatter` as a `__FN__` sentinel INSIDE the skeleton
   returned by `serialize()`, unlike every other type (whose skeletons never carry
   formatter functions, because the interactive JS renderer attaches its own). This is
   safe here because the row->title mapping never changes for hidden/renamed/highlighted
   samples (only for a change in which METRICS exist), and because the interactive JS
   path (Task 2.2, not this module) simply overwrites `axisLabel.formatter` with its own
   live function before `setOption()`; only the SSR path (`static_export.py`'s `__FN__`
   walker) actually executes the sentinel this module emits.

2. Pre-computed polygons: `serialize()` additionally stores, per dataset, a `"violins"`
   key (`{metric: {"poly": [[x, y], ...], "range": [lo, hi]}}`) computed from ALL samples
   (not toolbox-filtered). KDE is too expensive to redo in JS on every render; the
   interactive renderer uses this polygon directly when no samples are hidden, and
   recomputes with the JS port of `kde()` only when samples are actually hidden. The
   inner Q1-Q3 box + median line (POLISH.md #10b) is NOT cached here: a sort + linear
   interpolation is cheap enough that the JS side (`buildSeries()`) just recomputes it
   from the currently-visible values on every render, same as the beeswarm points.

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
from multiqc.plots.echarts.converter import convert_layout
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

_DEFAULT_STROKE = (
    "rgb(31,77,158)"  # equivalent to the old "#1f4d9e", kept as rgb(...) so _rgba() below applies uniformly
)
_DEFAULT_FILL = f"rgba(91,143,249,{FILL_ALPHA})"
_DEFAULT_BOX_FILL = f"rgba(31,77,158,{BOX_FILL_ALPHA})"
_SCATTER_COLOR = "rgba(30,50,80,0.85)"


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
    Per-metric row normalization range. Prefers the column's own configured axis range
    (`header.xaxis.range`, the exact `[dmin, dmax]` Plotly's per-row x-axis uses, see
    `Dataset.create` in `multiqc/plots/violin.py`), so a violin occupies the same
    fraction of its row that Plotly's does instead of always stretching to fill the full
    row width regardless of how narrow the real data spread is (POLISH.md #10c). Falls
    back to a padded observed-value range (BUILD_PLAN.md Task 2.1), matching Plotly's own
    autorange, when the column has no configured range.
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
        # instead of pinning it to the row's left edge, so the (x - lo) / (hi - lo)
        # normalization below never divides by zero.
        lo, hi = values[0] - 0.5, values[0] + 0.5
    return lo, hi


def _metric_title(header: ViolinColumn) -> str:
    if header.namespace:
        return f"{header.namespace}: {header.title}"
    return header.title


def _violin_polygon(values: List[float], lo: float, hi: float, row_idx: int) -> List[List[float]]:
    """
    The KDE polygon for one row, in normalized (0..1) x-space and integer-offset y-space:
    a closed shape formed by the density curve above the row's baseline followed by the
    same curve, mirrored, below it. Ported verbatim from `scripts/06_violin_final.py` lines
    88-98.
    """
    span = hi - lo
    if max(values) == min(values):
        # Degenerate (zero-variance) metric: kde()'s Silverman bandwidth collapses to its
        # 1e-9 floor here, which would otherwise blow up into a single needle-thin,
        # full-height spike (POLISH.md #10d). Draw a flat row instead, matching Plotly's
        # own near-invisible rendering for constant-value metrics; the median tick drawn
        # over it (see series()) still marks the value.
        x = (values[0] - lo) / span
        return [[x, row_idx], [x, row_idx]]

    xs = [lo + span * i / (N_BINS - 1) for i in range(N_BINS)]
    ys = kde(values, xs)
    ymax = max(ys) or 1.0
    top = [[(x - lo) / span, row_idx + (y / ymax) * ROW_HEIGHT] for x, y in zip(xs, ys)]
    bottom = [[(x - lo) / span, row_idx - (y / ymax) * ROW_HEIGHT] for x, y in zip(reversed(xs), reversed(ys))]
    return top + bottom


def violin_polygons(dataset: Dataset) -> Dict[str, Dict[str, Any]]:
    """
    `{metric: {"poly": [[x, y], ...], "range": [lo, hi]}}`, computed from ALL samples in
    the dataset (not toolbox-filtered): the VIOLIN EXCEPTION described in the module
    docstring. `range` is the padded `[lo, hi]` domain used to normalize `poly`'s x
    coordinates (i.e. the same pair `_violin_polygon` used), so a consumer can denormalize
    a polygon x back to a real value via `lo + x * (hi - lo)`.
    """
    result: Dict[str, Dict[str, Any]] = {}
    for row_idx, metric in enumerate(_visible_metrics(dataset)):
        header = dataset.header_by_metric[metric]
        values = _numeric_values(dataset, metric)
        lo, hi = _metric_range(header, values)
        poly = _violin_polygon(values, lo, hi, row_idx)
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
    q1x: float,
    medx: float,
    q3x: float,
    row_idx: int,
    fill: str,
    stroke: str,
    box_fill: str,
) -> Dict[str, Any]:
    """
    One `custom` series per row whose `renderItem` draws the KDE polygon computed in
    Python, plus an inner Q1-Q3 box and median line (POLISH.md #10b), matching Plotly's
    `box_visible`/`meanline` violin config. `q1x`/`medx`/`q3x` are already normalized to
    the same 0..1 x-space as `poly`. This is the ONE place a real function needs to reach
    the SSR bundle: the `__FN__` sentinel body below is plain MultiQC-generated JS source
    (a JSON-encoded point list plus fixed numeric/color literals), never user data.
    """
    body = (
        f"var poly = {json.dumps(poly)};"
        "var pts = poly.map(function(p) { return api.coord(p); });"
        f"var bTL = api.coord([{json.dumps(q1x)}, {json.dumps(row_idx - BOX_HALF_HEIGHT)}]);"
        f"var bBR = api.coord([{json.dumps(q3x)}, {json.dumps(row_idx + BOX_HALF_HEIGHT)}]);"
        f"var mTop = api.coord([{json.dumps(medx)}, {json.dumps(row_idx - MEDIAN_HALF_HEIGHT)}]);"
        f"var mBot = api.coord([{json.dumps(medx)}, {json.dumps(row_idx + MEDIAN_HALF_HEIGHT)}]);"
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
        "data": [0],
        "renderItem": {"__FN__": True, "body": body},
        "silent": True,
        "z": 1,
    }


def _annotation_render_series(row_idx: int, lo_x: float, hi_x: float, lo_obs: float, hi_obs: float) -> Dict[str, Any]:
    """
    One `custom` series per row drawing the observed min/max text labels at the violin's
    actual extent within the row (`lo_x`/`hi_x`, normalized 0..1): rows no longer always
    span the full row width (POLISH.md #10c), so labels must follow the drawn blob rather
    than always sitting at the row's physical ends.
    """
    lo_label = json.dumps(f"{lo_obs:g}")
    hi_label = json.dumps(f"{hi_obs:g}")
    body = (
        f"var L = api.coord([{json.dumps(lo_x)}, {row_idx}]); var R = api.coord([{json.dumps(hi_x)}, {row_idx}]);"
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
        "data": [0],
        "renderItem": {"__FN__": True, "body": body},
        "silent": True,
        "z": 3,
    }


def _scatter_series_for_metric(
    dataset: Dataset, metric: ColumnAnchor, row_idx: int, lo: float, hi: float
) -> Dict[str, Any]:
    """
    ONE beeswarm `scatter` series for this metric (not one shared series across all
    metrics, unlike the reference script): `jitter`/`jitterOverlap` gives ECharts 6's
    native beeswarm spread. Items follow the exact shape asked for in
    `multiqc-echarts-exploration/BUILD_PLAN.md` Task 2.1: `{"value": [normx, rowIdx],
    "name": sample}`; the real (denormalized) value and per-sample coloring are added by
    the interactive JS builder in Task 2.2, not here.
    """
    header = dataset.header_by_metric[metric]
    data: List[Dict[str, Any]] = []
    if header.show_points:
        if header.show_only_outliers:
            values_by_sample = dataset.scatter_value_by_sample_by_metric.get(metric, {})
        else:
            values_by_sample = dataset.violin_value_by_sample_by_metric.get(metric, {})
        span = hi - lo
        for sample, value in values_by_sample.items():
            if not isinstance(value, (int, float)):
                continue
            normx = (float(value) - lo) / span
            data.append({"value": [normx, row_idx], "name": str(sample)})

    return {
        "type": "scatter",
        "name": str(metric),
        "data": data,
        "symbolSize": 6,
        "jitter": 22,
        "jitterOverlap": False,
        "itemStyle": {"color": _SCATTER_COLOR, "borderColor": "#fff", "borderWidth": 0.5},
        "z": 2,
    }


def layout_option(plot: "Plot[Any, Any]", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one violin dataset, minus the KDE/annotation/scatter
    `series` (built by `series()` below for the SSR/get_option path). See the module
    docstring for why, unlike every other type, this skeleton DOES carry a `__FN__`
    formatter sentinel on `yAxis.axisLabel` (the VIOLIN EXCEPTION).
    """
    option = convert_layout(plot.layout, dataset.layout)

    metrics = _visible_metrics(dataset)
    n = len(metrics)
    row_titles = {i: _metric_title(dataset.header_by_metric[m]) for i, m in enumerate(metrics)}

    # Hidden value x-axis: violin x-position is normalized per-row to 0..1 (with a small
    # margin so polygons/points near the row's observed min/max are not clipped).
    option["xAxis"] = {"type": "value", "min": -0.05, "max": 1.05, "show": False}

    # VALUE-AXIS TRICK (module docstring point 1): a category yAxis snaps each violin's
    # fractional offset to the nearest integer tick, flattening every polygon. Using a
    # value axis with integer positions and a row-index -> metric-title formatter avoids
    # that, at the cost of needing a formatter function in the skeleton (see docstring).
    option["yAxis"] = {
        "type": "value",
        "inverse": True,
        "min": -0.5,
        "max": max(n - 0.5, 0.5),
        "interval": 1,
        "axisLabel": {
            "fontSize": 12,
            "align": "right",
            "verticalAlign": "middle",
            "formatter": {
                "__FN__": True,
                "body": f"var m = {json.dumps(row_titles)}; return m[Math.round(v)] || '';",
            },
        },
        "axisTick": {"show": False},
        "axisLine": {"show": False},
        "splitLine": {"show": False},
    }
    option["tooltip"]["trigger"] = "item"

    # No cross spike lines (POLISH.md #13): both axes here are the VALUE-AXIS TRICK
    # (module docstring), not real data coordinates (x is normalized 0..1 per row, y is
    # a continuous row index), so a value label following the cursor would show numbers
    # that mean nothing to the viewer. The crosshair mouse *cursor* (set globally in JS,
    # see echarts-plotting.js's buildCurrentOption) still applies to violin like every
    # other type; only the axis-tracking guide lines are skipped here.
    option["tooltip"]["axisPointer"] = {"show": False}
    option["axisPointer"] = {"show": False}

    return option


def series(dataset: Dataset, pconfig: TableConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    `n_metric` violin (KDE polygon + inner box + median) `custom` series, then `n_metric`
    min/max-annotation `custom` series, then `n_metric` beeswarm `scatter` series (one per
    metric; see `_scatter_series_for_metric`). This is the SSR/get_option (non-toolbox)
    path; the interactive path is `EchartsViolinPlot.buildSeries()`
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
        span = hi - lo

        poly = _violin_polygon(values, lo, hi, row_idx)
        fill, stroke, box_fill = _violin_colors(header)
        q1, median, q3 = _row_quartiles(values)
        violin_series.append(
            _violin_render_series(
                poly, (q1 - lo) / span, (median - lo) / span, (q3 - lo) / span, row_idx, fill, stroke, box_fill
            )
        )

        lo_obs, hi_obs = min(values), max(values)
        annotation_series.append(
            _annotation_render_series(row_idx, (lo_obs - lo) / span, (hi_obs - lo) / span, lo_obs, hi_obs)
        )

        scatter_series.append(_scatter_series_for_metric(dataset, metric, row_idx, lo, hi))

    return violin_series + annotation_series + scatter_series


def axis_data(dataset: Dataset, pconfig: TableConfig) -> Optional[List[Tuple[str, List[str]]]]:
    """
    Always `None`: both violin axes are `type: "value"` (the x-axis is a hidden
    normalized 0..1 axis, the y-axis is the VALUE-AXIS TRICK), so there is no category
    `data` array for the toolbox to fill in.
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
