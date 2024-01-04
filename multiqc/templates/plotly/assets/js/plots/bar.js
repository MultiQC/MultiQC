class BarPlot extends Plot {
  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    let cats = this.datasets[this.active_dataset_idx]["cats"];
    if (cats.length === 0) return 0; // no categories
    return cats[0].data.length; // no data for a category
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let cats = this.datasets[this.active_dataset_idx]["cats"];
    let samples = this.datasets[this.active_dataset_idx]["samples"];
    if (cats === 0 || samples === 0) return [];

    let samplesSettings = applyToolboxSettings(samples);
    if (samplesSettings == null) return []; // All series are hidden, do not render the graph.

    // Rename and filter samples:
    let filteredSettings = samplesSettings.filter((s) => !s.hidden);
    samples = filteredSettings.map((s) => s.name);

    return cats.map((cat) => {
      let data = this.p_active ? cat["data_pct"] : cat.data;
      data = data.filter((_, si) => !samplesSettings[si].hidden);

      let params = JSON.parse(JSON.stringify(this.trace_params)); // deep copy
      params.marker.color = cat.color;
      params.marker.line = {
        // Remove grey from highlights:
        color: filteredSettings.map((s) => (s.highlight && s.highlight !== "#cccccc" ? s.highlight : null)),
        width: filteredSettings.map((s) => (s.highlight && s.highlight !== "#cccccc" ? 2 : params.marker.line.width)),
      };

      return {
        type: "bar",
        x: data,
        y: samples,
        name: cat.name,
        ...params,
      };
    });
  }
}
