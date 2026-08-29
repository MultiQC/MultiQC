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
   recomputes with the JS port of `kde()` only when samples are actually hidden.

Trust boundary: every `{"__FN__": True, "body": "<js source>"}` sentinel emitted in this
module is built from typed, MultiQC-generated values (floats, sample names formatted with
`json.dumps`), never raw external input; `static_export.py` turns these bodies into real
JS functions via `new Function(...)`, which is only safe because of that provenance (see
the comment at `static_export._FN_WALKER_JS`).
"""

import json
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Tuple

from multiqc.plots.echarts.converter import convert_layout
from multiqc.plots.table_object import TableConfig
from multiqc.plots.violin import ColumnAnchor, Dataset, ViolinColumn
from multiqc.utils.mqc_colour import color_to_rgb_string

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot

ROW_HEIGHT = 0.42
N_BINS = 60
RANGE_PAD = 0.15

_DEFAULT_FILL = "rgba(91,143,249,0.5)"
_DEFAULT_STROKE = "#1f4d9e"
_SCATTER_COLOR = "rgba(30,50,80,0.85)"


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


def _metric_range(values: List[float]) -> Tuple[float, float]:
    """Per-metric row normalization range with 15% padding, per BUILD_PLAN.md Task 2.1."""
    lo, hi = min(values), max(values)
    span = hi - lo
    lo -= span * RANGE_PAD
    hi += span * RANGE_PAD
    if hi <= lo:
        # All values identical (span == 0): give the row a non-degenerate width so the
        # (x - lo) / (hi - lo) normalization below never divides by zero.
        hi = lo + 1.0
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
        values = _numeric_values(dataset, metric)
        lo, hi = _metric_range(values)
        poly = _violin_polygon(values, lo, hi, row_idx)
        result[str(metric)] = {"poly": poly, "range": [lo, hi]}
    return result


def _violin_colors(header: ViolinColumn) -> Tuple[str, str]:
    if not header.color:
        return _DEFAULT_FILL, _DEFAULT_STROKE
    stroke = color_to_rgb_string(header.color)  # "rgb(r,g,b)"
    fill = stroke.replace("rgb(", "rgba(").replace(")", ",0.5)")
    return fill, stroke


def _kde_render_series(poly: List[List[float]], fill: str, stroke: str) -> Dict[str, Any]:
    """
    One `custom` series per row whose `renderItem` draws the KDE polygon computed in
    Python. This is the ONE place a real function needs to reach the SSR bundle: the
    `__FN__` sentinel body below is plain MultiQC-generated JS source (a JSON-encoded
    point list plus fixed fill/stroke strings), never user data.
    """
    body = (
        f"var poly = {json.dumps(poly)};"
        "var pts = poly.map(function(p) { return api.coord(p); });"
        "return { type: 'polygon', shape: { points: pts }, "
        f"style: {{ fill: {json.dumps(fill)}, stroke: {json.dumps(stroke)}, lineWidth: 1 }} }};"
    )
    return {
        "type": "custom",
        "coordinateSystem": "cartesian2d",
        "data": [0],
        "renderItem": {"__FN__": True, "body": body},
        "silent": True,
        "z": 1,
    }


def _annotation_render_series(row_idx: int, lo_obs: float, hi_obs: float) -> Dict[str, Any]:
    """One `custom` series per row drawing the observed min/max text labels at row ends."""
    lo_label = json.dumps(f"{lo_obs:g}")
    hi_label = json.dumps(f"{hi_obs:g}")
    body = (
        f"var L = api.coord([0, {row_idx}]); var R = api.coord([1, {row_idx}]);"
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

    return option


def series(dataset: Dataset, pconfig: TableConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    `n_metric` KDE-polygon `custom` series, then `n_metric` min/max-annotation `custom`
    series, then `n_metric` beeswarm `scatter` series (one per metric; see
    `_scatter_series_for_metric`). This is the SSR/get_option (non-toolbox) path; the
    interactive path is `EchartsViolinPlot.buildSeries()`
    (`templates/echarts/src/js/plots/violin.js`, Task 2.2).

    `pconfig`/`is_pct` are accepted for dispatch-signature parity with `bar.series`; violin
    plots have no percentage switch, so `is_pct` is unused here.
    """
    metrics = _visible_metrics(dataset)

    kde_series = []
    annotation_series = []
    scatter_series = []
    for row_idx, metric in enumerate(metrics):
        header = dataset.header_by_metric[metric]
        values = _numeric_values(dataset, metric)
        lo, hi = _metric_range(values)

        poly = _violin_polygon(values, lo, hi, row_idx)
        fill, stroke = _violin_colors(header)
        kde_series.append(_kde_render_series(poly, fill, stroke))

        annotation_series.append(_annotation_render_series(row_idx, min(values), max(values)))

        scatter_series.append(_scatter_series_for_metric(dataset, metric, row_idx, lo, hi))

    return kde_series + annotation_series + scatter_series


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
