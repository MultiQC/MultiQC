"""
ECharts option builder for `PlotType.BOX` (`multiqc/plots/box.py`).

MultiQC box plots are horizontal: samples on the (inverted) category axis (y), values
on the value axis (x), mirroring bar.py. See `templates/echarts/src/js/plots/box.js` for
the interactive (toolbox-aware) JS counterpart of `series()`/`axis_data()` below.

ECharts' native `boxplot` series type takes a PRECOMPUTED five-number summary
(`[min, q1, median, q3, max]`) per box; unlike Plotly's `go.Box`, it never computes
quartiles from raw points itself. Outlier points (values outside the box whiskers) are
drawn as a companion `scatter` series, the same convention ECharts' own
`echarts.dataTool.prepareBoxplotData` helper uses.

Quartile method (risk flagged in `multiqc-echarts-exploration/BUILD_PLAN.md` section 4,
"Box quartile formulas"): linear interpolation between the two closest ranks, i.e.
numpy's default `np.percentile` method and R's type-7 quantile, which is also Plotly's
default `quartilemethod="linear"` for `go.Box`. Whisker ends are the most extreme
in-fence data points (Tukey convention, 1.5 * IQR), matching both Plotly's
`lowerfence`/`upperfence` and ECharts' own boxplot data-prep tool; so `min`/`max` in the
five-number summary below are NOT the sample's absolute min/max, but the whisker ends;
values outside the fence are outliers.

This is the SAME formula ported verbatim to JS in
`templates/echarts/src/js/plots/box.js` (`quantile`/`fiveNumberSummary`/`outliers`); the
cross-language contract is pinned by one golden fixed-input/expected-output test in
`tests/test_plots_echarts.py` (`test_box_five_number_summary_golden_values` /
`test_box_outliers_golden_values`), whose exact input and expected values are mirrored
in a comment block at the top of that JS file.

For `is_stats_data` datasets (`Dataset` in `multiqc/plots/box.py`, pre-calculated
per-sample statistics, used for memory efficiency on huge datasets), the five-number
summary is read directly from the sample's `min`/`q1`/`median`/`q3`/`max` keys: we do
not (and, lacking raw values, cannot) recompute quartiles or outliers for these samples.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Tuple, Union, cast

from multiqc.plots.box import BoxPlotConfig, BoxStatsT, Dataset
from multiqc.plots.echarts.converter import convert_layout, zoom_option
from multiqc.utils.mqc_colour import color_to_rgb_string

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot

_DEFAULT_COLOR = "#4899e8"  # matches multiqc/plots/box.py's trace_params marker color


def _quantile(sorted_values: List[float], q: float) -> float:
    """
    Linear-interpolation quantile (numpy's default `np.percentile` method, R's type-7,
    and Plotly's default `quartilemethod="linear"`). `sorted_values` must already be
    sorted ascending and non-empty.
    """
    n = len(sorted_values)
    if n == 1:
        return sorted_values[0]
    h = (n - 1) * q
    lo = int(h)  # floor: h is always >= 0
    hi = min(lo + 1, n - 1)
    frac = h - lo
    return sorted_values[lo] + frac * (sorted_values[hi] - sorted_values[lo])


def _fences(q1: float, q3: float) -> Tuple[float, float]:
    iqr = q3 - q1
    return q1 - 1.5 * iqr, q3 + 1.5 * iqr


def five_number_summary(values: List[Union[int, float]]) -> List[float]:
    """
    `[min, q1, median, q3, max]` for one sample's raw data points, in the format the
    ECharts `boxplot` series expects. `min`/`max` here are the boxplot WHISKER ends (the
    most extreme data points within 1.5 * IQR of the box, Tukey convention), not the
    sample's absolute min/max; see the module docstring.
    """
    sorted_values = sorted(float(v) for v in values)
    q1 = _quantile(sorted_values, 0.25)
    median = _quantile(sorted_values, 0.5)
    q3 = _quantile(sorted_values, 0.75)
    lower_fence, upper_fence = _fences(q1, q3)
    in_fence = [v for v in sorted_values if lower_fence <= v <= upper_fence]
    whisker_min = in_fence[0] if in_fence else sorted_values[0]
    whisker_max = in_fence[-1] if in_fence else sorted_values[-1]
    return [whisker_min, q1, median, q3, whisker_max]


def outliers(values: List[Union[int, float]]) -> List[float]:
    """Values outside the Tukey 1.5 * IQR fence: the points the whisker box doesn't reach."""
    sorted_values = sorted(float(v) for v in values)
    q1 = _quantile(sorted_values, 0.25)
    q3 = _quantile(sorted_values, 0.75)
    lower_fence, upper_fence = _fences(q1, q3)
    return [v for v in sorted_values if v < lower_fence or v > upper_fence]


def _stats_five_number(stats: BoxStatsT) -> List[float]:
    return [float(stats["min"]), float(stats["q1"]), float(stats["median"]), float(stats["q3"]), float(stats["max"])]


def _box_item_style(color: str) -> Dict[str, str]:
    """
    Match Plotly's own box rendering exactly (verified against the default template's
    output: a single stroked SVG path per box at `stroke-opacity: 1` for the border,
    whiskers AND median line, with `fill-opacity: 0.5` for the body). ECharts draws the
    median line using `itemStyle.borderColor`, so a solid border against a lighter,
    semi-transparent fill is what makes the median visible instead of blending into a
    fully opaque box.
    """
    border = color_to_rgb_string(color)  # "rgb(r,g,b)"
    fill = border.replace("rgb(", "rgba(").replace(")", ",0.5)")
    return {"color": fill, "borderColor": border}


def layout_option(plot: "Plot[Any, Any]", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one box-plot dataset, minus `series` and axis
    `data` arrays (samples are toolbox-dependent; the JS renderer fills `yAxis.data`).
    """
    # xAxis is the value axis (box plots are horizontal); Plotly-style data-fitted
    # autorange instead of ECharts' forced-0 default, see `converter._convert_axis`.
    option = convert_layout(plot.layout, dataset.layout, scale_x=True)

    option["yAxis"]["type"] = "category"
    option["yAxis"]["inverse"] = True
    option["tooltip"]["trigger"] = "item"
    # Plotly-style click+drag box-zoom (POLISH.md #17), X (value) axis only: the y-axis
    # is a sample-name category list already managed by MultiQC's own sidebar toolbox
    # (hide/highlight), see `zoom_option`.
    option.update(zoom_option(x=True, y=False))

    return option


def _build(dataset: Dataset) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Shared computation for `series()` and `mark_count()`, so the two can never disagree
    on how many marks are actually drawn: the boxplot five-number data (one entry per
    sample, in `dataset.samples` order) and the companion outlier/point scatter data.
    """
    boxpoints = dataset.trace_params.get("boxpoints", "outliers")
    base_color = dataset.trace_params.get("marker", {}).get("color") or _DEFAULT_COLOR
    item_style = _box_item_style(base_color)

    box_data: List[Dict[str, Any]] = []
    scatter_data: List[Dict[str, Any]] = []
    for yi, (sample, values) in enumerate(zip(dataset.samples, dataset.data)):
        if dataset.is_stats_data:
            box_data.append({"value": _stats_five_number(cast(BoxStatsT, values)), "itemStyle": item_style})
            continue

        raw_values = cast(List[Union[int, float]], values)
        if not raw_values:
            box_data.append({"value": [0.0, 0.0, 0.0, 0.0, 0.0], "itemStyle": item_style})
            continue
        box_data.append({"value": five_number_summary(raw_values), "itemStyle": item_style})

        if boxpoints is False:
            continue
        shown = sorted(float(v) for v in raw_values) if boxpoints == "all" else outliers(raw_values)
        scatter_data.extend({"name": sample, "value": [value, yi]} for value in shown)

    return box_data, scatter_data


def series(dataset: Dataset, pconfig: BoxPlotConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    One `{"type": "boxplot"}` series (one five-number entry per sample, positionally
    matching `yAxis.data`) plus, unless `boxpoints` is `False`, a companion
    `{"type": "scatter"}` series of outlier points (or every point, when
    `boxpoints == "all"`). This is the SSR/get_option (non-toolbox) path; the
    interactive path is `EchartsBoxPlot.buildSeries()`
    (`templates/echarts/src/js/plots/box.js`).

    `pconfig`/`is_pct` are accepted for dispatch-signature parity with `bar.series`; box
    plots never enable the percentage switch, so they are unused here.
    """
    box_data, scatter_data = _build(dataset)
    result: List[Dict[str, Any]] = [{"type": "boxplot", "data": box_data}]
    if scatter_data:
        result.append({"type": "scatter", "name": "outliers", "data": scatter_data, "symbolSize": 5})
    return result


def axis_data(dataset: Dataset, pconfig: BoxPlotConfig) -> List[Tuple[str, List[str]]]:
    """`[("yAxis", sample names)]`: box plots are always horizontal, samples on `yAxis`."""
    return [("yAxis", list(dataset.samples))]


def mark_count(dataset: Dataset) -> int:
    """One mark per sample (box), plus one per outlier/shown point."""
    box_data, scatter_data = _build(dataset)
    return len(box_data) + len(scatter_data)
