class BarPlot extends Plot {
  activeDatasetSamples() {
    if (this.datasets.length === 0) return [];
    let ds = this.datasets[this.active_dataset_idx];
    if (ds.length === 0) return []; // no categories
    let cat = ds[0]; // samples should be same in every category
    return cat["samples"];
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let dataset = this.datasets[this.active_dataset_idx];

    let samples = applyToolboxSettings(this.activeDatasetSamples());

    // All series are hidden, do not render the graph.
    if (samples == null) return;

    // Rename samples
    dataset.each((cat) => cat["samples"].each((sn, i) => (cat["samples"][i] = samples[sn].name)));

    // Filter out hidden samples
    let visibleSamples = samples.filter((sn) => !samples[sn].hidden);

    return dataset.map((cat) => {
      let data = this.p_active ? cat.data_pct : cat.data;
      let visibleData = data.filter((x, i) => !samples[d["samples"][i]].hidden);
      return {
        type: "bar",
        x: visibleData,
        y: visibleSamples,
        name: cat.name,
        orientation: "h",
        marker: {
          line: {
            // Remove grey from highlights, as we don't need to remove default coloring of background samples
            color: samples.map((x) => (x.highlight && x.highlight !== "#cccccc" ? x.highlight : null)),
            width: samples.map((x) => (x.highlight && x.highlight !== "#cccccc" ? 2 : 0)),
          },
        },
        marker_color: cat.color,
      };
    });
  }

  // TODO: perhaps use it or Plotly.react
  // replot() {
  //   // Updates plot given new underlying data
  //   let dataset = this.datasets[this.active_dataset_idx];
  //   let x = dataset.map(
  //     (cat) => this.p_active ? cat.data_pct : cat.data
  //   );
  //   Plotly.restyle(this.target, "x", x);
  // };
}
