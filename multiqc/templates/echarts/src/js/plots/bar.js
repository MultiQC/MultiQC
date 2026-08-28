// ECharts bar plot: extends the default template's BarPlot (imported just before this
// file in main-js.js) purely for its prepData()/exportData()/formatDatasetForAiPrompt()/
// activeDatasetSize(), which are engine-neutral (they only read/filter raw dataset
// fields through applyToolboxSettings()). window.BarPlot itself extends window.Plot,
// which by this point in the import order is OUR echarts Plot class (from
// echarts-plotting.js), not Plotly's.
class EchartsBarPlot extends window.BarPlot {
  // Constructs and returns ECharts bar series for the active dataset, applying the
  // toolbox (hide/rename/highlight) and pct toggle via the inherited prepData().
  buildSeries() {
    let dataset = this.datasets[this.activeDatasetIdx];

    // Multicategory (sample groups) is out of scope for Phase 0 (see BUILD_PLAN.md
    // Task 0.5); default's prepData() would silently take its Plotly-only multicategory
    // branch instead, so we crash loudly here rather than let that happen.
    if (dataset["group_labels"]) {
      throw new Error("ECharts bar plot does not support sample groups yet");
    }

    let [cats] = this.prepData();

    if (cats.length === 0 || this.filteredSettings.length === 0) {
      this._axisData = [];
      return [];
    }

    // Stashed here for renderPlot() to assign to option.yAxis.data (samples are
    // toolbox-dependent, so the skeleton from Python never contains them).
    this._axisData = this.filteredSettings.map((s) => s.name);

    let highlighted = this.filteredSettings.filter((s) => s.highlight);
    let barmode = this.echarts.datasets[this.activeDatasetIdx].layout["_mqc"]?.barmode;
    let isGroup = barmode === "group";

    return cats.map((cat) => {
      let series = {
        type: "bar",
        name: cat.name,
        barCategoryGap: "30%",
        data: this.filteredSettings.map((sample, sampleIdx) => {
          // Dim non-highlighted samples to 0.1 alpha when any highlight is active,
          // same rule as default plots/bar.js.
          let alpha = highlighted.length > 0 && sample.highlight === null ? 0.1 : 1;
          return {
            value: cat.data[sampleIdx],
            itemStyle: { color: "rgba(" + cat.color + "," + alpha + ")" },
          };
        }),
      };
      if (!isGroup) series.stack = "total";
      return series;
    });
  }
}

// Make EchartsBarPlot globally available (referenced by bare window.EchartsBarPlot in
// echarts-plotting.js's initPlot(), not a bare identifier, so it's robust across the
// minified bundle regardless of module ordering after bundling).
window.EchartsBarPlot = EchartsBarPlot;
