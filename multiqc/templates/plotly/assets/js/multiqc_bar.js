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
    let samples = this.activeDatasetSamples();

    // Rename samples
    renameSamples(samples, function (sample_idx, new_name) {
      for (let cat of dataset) cat["samples"][sample_idx] = new_name;
    });

    // Hide samples
    let plot_group_div = $("#" + this.target).closest(".mqc_hcplot_plotgroup");
    let idx_to_hide = hideSamples(plot_group_div, samples);
    if (idx_to_hide.length === samples.length)
      // All series are hidden. Do not render the graph.
      return;

    // Filter out hidden samples
    let visible_samples = samples.filter((x, i) => !idx_to_hide.includes(i));

    // Highlight samples
    let highlight_colors = getHighlightColors(visible_samples);
    // Remove grey, as we don't need to remove default coloring of background samples
    highlight_colors = highlight_colors.map((x) => (x === "#cccccc" ? null : x));

    // Render the plotly plot
    let traces = [];
    for (let cat of dataset) {
      let data = this.p_active ? cat.data_pct : cat.data;
      let visible_data = data.filter((x, i) => !idx_to_hide.includes(i));
      let trace = {
        type: "bar",
        x: visible_data,
        y: visible_samples,
        name: cat.name,
        orientation: "h",
        marker: {
          line: {
            color: highlight_colors,
            width: highlight_colors.map((x) => (x ? 2 : 0)),
          },
        },
        marker_color: cat.color,
      };
      traces.push(trace);
    }
    return traces;
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
