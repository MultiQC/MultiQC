// ECharts scatter plot: extends the default template's ScatterPlot (imported just before
// this file in main-js.js) purely for its prepData()/exportData()/formatDatasetForAiPrompt()/
// activeDatasetSize(), which are engine-neutral (they only read/filter raw dataset fields
// through applyToolboxSettings()). window.ScatterPlot itself extends window.Plot, which by
// this point in the import order is OUR echarts Plot class (from echarts-plotting.js), not
// Plotly's.
//
// Series strategy: ONE `{type: "scatter"}` series per dataset with per-item styling, not
// one series per point or per color-group (see multiqc/plots/echarts/scatter.py docstring
// for why: it keeps a large point cloud a single series, which is what makes the
// canvas-vs-svg renderer threshold meaningful).
class EchartsScatterPlot extends window.ScatterPlot {
  buildSeries() {
    let dataset = this.datasets[this.activeDatasetIdx];
    let marker = dataset["trace_params"]?.["marker"] || {};
    let defaultColor = marker.color ?? "rgba(124, 181, 236, .5)";
    let defaultSize = marker.size ?? 10;

    let [samples, points] = this.prepData();
    // prepData() (inherited from default ScatterPlot) leaves a hole (undefined) for
    // every point the toolbox hides instead of removing it; drop those here.
    points = points.filter(Boolean);

    if (points.length === 0 || samples.length === 0) return [];

    // Reorder so highlighted points draw on top, same as default plots/scatter.js.
    let highlighted = points.filter((p) => p.highlight);
    let nonHighlighted = points.filter((p) => !p.highlight);
    points = nonHighlighted.concat(highlighted);

    let data = points.map((point) => {
      // Same highlight rule as default plots/scatter.js: once any point is
      // highlighted, recolor every point to its highlight color (or gray), replacing
      // whatever per-point color it had.
      let color = highlighted.length > 0 ? (point.highlight ?? "#cccccc") : (point.color ?? defaultColor);
      let itemStyle = { color: color };
      if (point.opacity !== undefined) itemStyle.opacity = point.opacity;
      if (point.marker_line_width !== undefined) itemStyle.borderWidth = point.marker_line_width;

      return {
        value: [point.x, point.y],
        name: point.name,
        symbolSize: point.marker_size ?? defaultSize,
        itemStyle: itemStyle,
      };
    });

    return [
      {
        type: "scatter",
        name: dataset["label"] ?? "",
        data: data,
      },
    ];
  }

  // Tooltip formatter can't live in the JSON-safe serialized skeleton, so it's attached
  // here to the fully-assembled option (see Plot.applyOptionOverrides in
  // echarts-plotting.js).
  applyOptionOverrides(option) {
    if (!option.tooltip) return;
    let xsuffix = this.layout?.xaxis?.ticksuffix ?? "";
    let ysuffix = this.layout?.yaxis?.ticksuffix ?? "";
    option.tooltip.formatter = (params) => {
      let [x, y] = params.value;
      return `<b>${params.name}</b><br/>X: ${x}${xsuffix}<br/>Y: ${y}${ysuffix}`;
    };
  }
}

// Make EchartsScatterPlot globally available (referenced by bare window.EchartsScatterPlot
// in echarts-plotting.js's initPlot(), not a bare identifier, so it's robust across the
// minified bundle regardless of module ordering after bundling).
window.EchartsScatterPlot = EchartsScatterPlot;
