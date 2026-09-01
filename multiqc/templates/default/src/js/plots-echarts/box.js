// ECharts box plot: extends the default template's BoxPlot (imported just before this
// file in main-js.js) for its constructor/prepData()/exportData()/
// formatDatasetForAiPrompt()/activeDatasetSize(), which are engine-neutral: prepData()
// already resolves the sort-switch (dataset.data_sorted/samples_sorted) and applies
// toolbox hide/rename via applyToolboxSettings(). window.BoxPlot itself extends
// window.Plot, which by this point in the import order is OUR echarts Plot class (from
// echarts-plotting.js), not Plotly's. The default box.js sort-toggle click handler at
// its bottom (`renderPlot(anchor)`) is engine-neutral and reused as-is.
//
// GOLDEN quartile fixture (cross-language contract with
// multiqc/plots/echarts/box.py, asserted in
// tests/test_plots_echarts.py::test_box_five_number_summary_golden_values /
// test_box_outliers_golden_values):
//   input:            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100]
//   fiveNumberSummary: [1, 3.5, 6, 8.5, 10]
//   outliers:          [100]
// fiveNumberSummary()/outliers() below MUST stay byte-for-byte the same formula as the
// Python functions of the same name in multiqc/plots/echarts/box.py.

// Linear-interpolation quantile (numpy's default `np.percentile` method, R's type-7,
// and Plotly's default `quartilemethod="linear"`); mirrors
// multiqc/plots/echarts/box.py::_quantile. `sortedValues` must already be sorted
// ascending and non-empty.
function quantile(sortedValues, q) {
  const n = sortedValues.length;
  if (n === 1) return sortedValues[0];
  const h = (n - 1) * q;
  const lo = Math.floor(h);
  const hi = Math.min(lo + 1, n - 1);
  const frac = h - lo;
  return sortedValues[lo] + frac * (sortedValues[hi] - sortedValues[lo]);
}

// [min, q1, median, q3, max] as ECharts' boxplot series expects, mirroring
// multiqc/plots/echarts/box.py::five_number_summary. min/max here are the boxplot
// WHISKER ends (most extreme in-fence points, Tukey 1.5*IQR convention), not the
// sample's absolute min/max; values outside are outliers (see outliers() below).
function fiveNumberSummary(values) {
  let sorted = [...values].sort((a, b) => a - b);
  let q1 = quantile(sorted, 0.25);
  let median = quantile(sorted, 0.5);
  let q3 = quantile(sorted, 0.75);
  let iqr = q3 - q1;
  let lowerFence = q1 - 1.5 * iqr;
  let upperFence = q3 + 1.5 * iqr;
  let inFence = sorted.filter((v) => v >= lowerFence && v <= upperFence);
  let whiskerMin = inFence.length ? inFence[0] : sorted[0];
  let whiskerMax = inFence.length ? inFence[inFence.length - 1] : sorted[sorted.length - 1];
  return [whiskerMin, q1, median, q3, whiskerMax];
}

// Values outside the Tukey 1.5*IQR fence, mirroring multiqc/plots/echarts/box.py::outliers.
function outliers(values) {
  let sorted = [...values].sort((a, b) => a - b);
  let q1 = quantile(sorted, 0.25);
  let q3 = quantile(sorted, 0.75);
  let iqr = q3 - q1;
  let lowerFence = q1 - 1.5 * iqr;
  let upperFence = q3 + 1.5 * iqr;
  return sorted.filter((v) => v < lowerFence || v > upperFence);
}

// [r, g, b] from any CSS color (hex, named like "grey", or an existing rgb()/rgba()
// string). Hex/named colors are resolved via a throwaway canvas context, which always
// normalizes an opaque color to "#rrggbb".
function colorToRgbTriplet(color) {
  let match = color.match(/rgba?\(([^)]+)\)/);
  if (match)
    return match[1]
      .split(",")
      .slice(0, 3)
      .map((s) => s.trim());
  let hex = Object.assign(document.createElement("canvas").getContext("2d"), { fillStyle: color }).fillStyle;
  return [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16));
}

// Match Plotly's own box rendering exactly (verified against the default template's
// output: a single stroked SVG path per box at stroke-opacity 1 for the border,
// whiskers AND median line, with fill-opacity 0.5 for the body). ECharts draws the
// median line using itemStyle.borderColor, so a solid border against a lighter,
// semi-transparent fill is what makes the median visible instead of blending into a
// fully opaque box. Mirrors multiqc/plots/echarts/box.py::_box_item_style.
function boxItemStyle(color) {
  let [r, g, b] = colorToRgbTriplet(color);
  return { color: `rgba(${r},${g},${b},0.5)`, borderColor: `rgb(${r},${g},${b})` };
}

class EchartsBoxPlot extends window.BoxPlot {
  // Constructs and returns the ECharts boxplot series (+ companion outlier scatter
  // series) for the active dataset, applying the sort switch and toolbox
  // (hide/rename/highlight) via the inherited prepData().
  buildSeries() {
    let [data, samples] = this.prepData();
    if (data.length === 0 || samples.length === 0) return [];

    let dataset = this.datasets[this.activeDatasetIdx];
    let traceParams = dataset["trace_params"] || {};
    let boxpoints = "boxpoints" in traceParams ? traceParams["boxpoints"] : "outliers";
    let baseColor = (traceParams["marker"] && traceParams["marker"]["color"]) || "#4899e8";

    // axisLabel (item B): colors the sample-name axis tick labels, mirroring Plotly's
    // recalculateTicks(); see Plot.sampleAxisLabel() in echarts-plotting.js. maxTicks
    // mirrors the Plotly box plot's own formula (templates/plotly/src/js/plots/box.js).
    let maxTicks = (this.layout.height - 140) / 12;
    this._axisData = {
      axis: "yAxis",
      data: this.filteredSettings.map((s) => s.name),
      axisLabel: this.sampleAxisLabel(this.filteredSettings, maxTicks),
    };

    let boxData = [];
    let scatterData = [];

    this.filteredSettings.forEach((sample, i) => {
      // Same highlight rule as default box.js's Plotly buildTraces(): recolor the
      // highlighted sample to its highlight color, every other sample to grey, once
      // ANY highlight is active anywhere (item A: global flag, not the local
      // `highlighted.length` count, so this plot greys out fully even with zero matches).
      let color = window.mqc_highlight_f_texts.length > 0 ? (sample.highlight ?? "grey") : baseColor;
      let values = data[i];

      if (this.isStatsData) {
        boxData.push({
          value: [values.min ?? 0, values.q1 ?? 0, values.median ?? 0, values.q3 ?? 0, values.max ?? 0],
          itemStyle: boxItemStyle(color),
        });
        return;
      }

      boxData.push({ value: fiveNumberSummary(values), itemStyle: boxItemStyle(color) });

      if (boxpoints === false) return;
      let shown = boxpoints === "all" ? [...values].sort((a, b) => a - b) : outliers(values);
      shown.forEach((value) => {
        scatterData.push({ name: sample.name, value: [value, i], itemStyle: { color: color } });
      });
    });

    let series = [{ type: "boxplot", layout: "horizontal", data: boxData }];
    if (scatterData.length > 0) {
      series.push({ type: "scatter", name: "outliers", data: scatterData, symbolSize: 5 });
    }
    return series;
  }

  // Tooltip formatter can't live in the JSON-safe serialized skeleton, so it's attached
  // here to the fully-assembled option (see Plot.applyOptionOverrides in
  // echarts-plotting.js). Shows the full five-number summary for a box, or the single
  // value for an outlier point.
  applyOptionOverrides(option) {
    if (!option.tooltip) return;
    let suffix = this.layout?.xaxis?.ticksuffix ?? "";
    // Box's value axis is xAxis (boxes are horizontal, samples on the category yAxis),
    // matching the Plotly reference (bar.js's own `layout.xaxis?.hoverformat` read applies
    // equally here since box.py also lays out horizontally). NOTE: read off the active
    // dataset's raw layout dict, not `this.layout` (see plots-echarts/bar.js's matching
    // comment): the neutral `AxisIR` behind `this.layout` drops `hoverformat`, a
    // Plotly-only key.
    let hoverformat = this.datasets[this.activeDatasetIdx]?.layout?.xaxis?.hoverformat;
    option.tooltip.formatter = (params) => {
      if (params.seriesType === "boxplot") {
        // ECharts prepends the category axis value/index to `value` for a boxplot on a
        // category axis (6 elements: [axisValue, min, q1, median, q3, max]), unlike the
        // plain 5-element five-number array in option.series[0].data; drop it here.
        let [minV, q1, median, q3, maxV] = params.value
          .slice(-5)
          .map((v) => window.formatWithHoverformat(v, hoverformat));
        return (
          `<b>${params.name}</b><br/>` +
          `Max: ${maxV}${suffix}<br/>Q3: ${q3}${suffix}<br/>Median: ${median}${suffix}<br/>` +
          `Q1: ${q1}${suffix}<br/>Min: ${minV}${suffix}`
        );
      }
      let value = Array.isArray(params.value) ? params.value[0] : params.value;
      return `<b>${params.name}</b>: ${window.formatWithHoverformat(value, hoverformat)}${suffix}`;
    };
  }
}

// Make EchartsBoxPlot globally available (referenced by bare window.EchartsBoxPlot in
// echarts-plotting.js's initPlot(), not a bare identifier, so it's robust across the
// minified bundle regardless of module ordering after bundling).
window.EchartsBoxPlot = EchartsBoxPlot;
