class BarPlot extends Plot {
  constructor(target, data) {
    super(target, data);
    for (let i = 0; i < data.datasets.length; i++) {
      data.datasets[i].samples = data.samples[i];
    }
  }

  active_dataset_size() {
    if (this.datasets.length === 0) return 0;
    return this.datasets[this.active_dataset_idx].samples.length;
  }

  // Constructs and returns traces for the Plotly plot
  build_traces() {
    let dataset = this.datasets[this.active_dataset_idx];

    // Rename samples
    rename_samples(dataset.samples, function (sample_idx, new_name) {
      dataset.samples[sample_idx] = new_name;
    });

    // Hide samples
    let plot_group_div = $("#" + this.target).closest(".mqc_hcplot_plotgroup");
    let idx_to_hide = hide_samples(plot_group_div, dataset.samples);
    if (idx_to_hide.length === dataset.samples.length)
      // All series are hidden. Do not render the graph.
      return;

    // Filter out hidden samples
    let visible_samples = dataset.samples.filter((x, i) => !idx_to_hide.includes(i));

    // Highlight samples
    let highlight_colors = get_highlight_colors(visible_samples);
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

  // replot() {
  //   // Updates plot given new underlying data
  //   let dataset = this.datasets[this.active_dataset_idx];
  //   let x = dataset.map(
  //     (cat) => this.p_active ? cat.data_pct : cat.data
  //   );
  //   Plotly.restyle(this.target, "x", x);
  // };
}
