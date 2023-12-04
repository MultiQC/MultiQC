class BarPlot extends Plot {
  constructor(data) {
    super(data);

    this.layout.showlegend = true;
  }

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

    let sampleNames = this.activeDatasetSamples();
    let samples = applyToolboxSettings(sampleNames);

    // All series are hidden, do not render the graph.
    if (samples == null) return;

    // Rename samples
    dataset.map((cat) => cat["samples"].map((sn, si) => (cat["samples"][si] = samples[si].name)));

    // Filter out hidden samples
    let visibleSamples = samples.filter((s) => !s.hidden);

    return dataset.map((cat) => {
      let data = this.p_active ? cat.data_pct : cat.data;
      let visibleData = data.filter((_, si) => !samples[si].hidden);
      return {
        type: "bar",
        x: visibleData,
        y: visibleSamples.map((s) => s.name),
        name: cat.name,
        orientation: "h",
        marker: {
          color: cat.color,
          line: {
            // Remove grey from highlights, as we don't need to remove default coloring of background samples
            color: visibleSamples.map((x) => (x.highlight && x.highlight !== "#cccccc" ? x.highlight : null)),
            width: visibleSamples.map((x) => (x.highlight && x.highlight !== "#cccccc" ? 2 : 0)),
          },
        },
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
