// ECharts line plot: extends the default template's LinePlot (imported just before this
// file in main-js.js) purely for its prepData()/exportData()/formatDatasetForAiPrompt()/
// activeDatasetSize(), which are engine-neutral (they only read/filter raw dataset
// fields through applyToolboxSettings()). window.LinePlot itself extends window.Plot,
// which by this point in the import order is OUR echarts Plot class (from
// echarts-plotting.js), not Plotly's.
class EchartsLinePlot extends window.LinePlot {
  // Constructs and returns ECharts line series for the active dataset, applying the
  // toolbox (hide/highlight) via the inherited prepData(), plus the static bands/lines
  // marker series (if any) stashed in the skeleton by multiqc/plots/echarts/line.py.
  buildSeries() {
    let dataset = this.datasets[this.activeDatasetIdx];
    let categorical = dataset["layout"]?.["xaxis"]?.["type"] === "category";
    // Whether markers are shown by default is a dataset-level decision (few enough total
    // points across all lines: linegraph.Dataset.create sets trace_params.mode), not a
    // per-series one; an explicit per-series line.marker (rare: e.g. extra_series
    // annotations) always shows its symbol regardless. Mirrors line.py::series().
    let mode = dataset["trace_params"]?.["mode"] || "";

    let [samples, lines] = this.prepData();

    if (lines.length === 0 || samples.length === 0) return [];

    // Reorder so highlighted lines draw on top, same as default plots/line.js.
    let highlighted = lines.filter((line) => line.highlight);
    let nonHighlighted = lines.filter((line) => !line.highlight);
    lines = nonHighlighted.concat(highlighted);

    let series = lines.map((line) => {
      // Same highlight rule as default plots/line.js: recolor non-highlighted lines to
      // gray once any line is highlighted, rather than bar's alpha-dimming (lines
      // overlap far more than stacked bars, so a flat gray reads better than alpha).
      let color = highlighted.length > 0 ? (line.highlight ?? "#cccccc") : line.color;
      return {
        type: "line",
        name: line.name,
        data: categorical ? line.pairs.map((p) => p[1]) : line.pairs.map((p) => [p[0], p[1]]),
        showSymbol: mode.includes("markers") || !!line.marker,
        smooth: false,
        lineStyle: { width: line.width, type: echartsDashType(line.dash) },
        itemStyle: { color: color },
        color: color,
      };
    });

    // Static quality-zone bands/lines (e.g. fastqc's green/yellow/red per-base-quality
    // background): precomputed once in Python, stashed in the skeleton, and here just
    // carried over onto a silent series so it survives the option.series = [...] reset.
    let bandsLines = this.echarts.datasets[this.activeDatasetIdx].layout["_mqc"]?.bandsLines;
    if (bandsLines) {
      series.push({
        type: "line",
        name: "",
        data: [],
        silent: true,
        showSymbol: false,
        tooltip: { show: false },
        ...bandsLines,
      });
    }

    return series;
  }

  // Tooltip formatter can't live in the JSON-safe serialized skeleton, so it's attached
  // here to the fully-assembled option (see Plot.applyOptionOverrides in
  // echarts-plotting.js).
  applyOptionOverrides(option) {
    if (!option.tooltip) return;
    option.tooltip.formatter = (params) => {
      let list = Array.isArray(params) ? params : [params];
      if (list.length === 0) return "";
      let xLabel = list[0].axisValueLabel ?? list[0].axisValue ?? list[0].name;
      let rows = list
        .map((p) => {
          let y = Array.isArray(p.value) ? p.value[1] : p.value;
          return `${p.marker}${p.seriesName}: <b>${window.formatNumber(y)}</b>`;
        })
        .join("<br/>");
      return `${xLabel}<br/>${rows}`;
    };
  }
}

// Plotly dash "dash"/"dot"/"dashdot"/"longdash"/"longdashdot"/"solid" (Series.dash is
// already normalized by convert_dash_style, multiqc/plots/plot.py) -> ECharts
// lineStyle.type, which only knows "solid"/"dashed"/"dotted". Mirrors the Python-side
// mapping in multiqc/plots/echarts/line.py::_echarts_dash (kept in sync manually: this
// is JS-only data, not shipped from Python for interactive lines built in JS).
function echartsDashType(dash) {
  const map = {
    solid: "solid",
    dash: "dashed",
    longdash: "dashed",
    dot: "dotted",
    dashdot: "dashed",
    longdashdot: "dashed",
  };
  return map[dash] ?? "solid";
}

// Make EchartsLinePlot globally available (referenced by bare window.EchartsLinePlot in
// echarts-plotting.js's initPlot(), not a bare identifier, so it's robust across the
// minified bundle regardless of module ordering after bundling).
window.EchartsLinePlot = EchartsLinePlot;
