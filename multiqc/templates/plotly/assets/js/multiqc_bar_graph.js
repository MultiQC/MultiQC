// Stacked Bar Graph
function plot_stacked_bar_graph(plot, target, dataset_idx) {
  if (plot === undefined || plot["plot_type"] !== "bar_graph") return false;

  if (dataset_idx === undefined) dataset_idx = 0;

  // Make a clone of everything, so that we can mess with it,
  // while keeping the original data intact
  let dataset = JSON.parse(JSON.stringify(plot["datasets"][dataset_idx]));
  let samples = JSON.parse(JSON.stringify(plot["samples"][dataset_idx]));
  let layout = JSON.parse(JSON.stringify(plot["layout"]));

  // Rename samples
  rename_samples(samples, function (sample_idx, new_name) {
    samples[sample_idx] = new_name;
  });

  // Highlight samples
  let highlight_colors = get_highlight_colors(samples);
  // Remove grey, as we don't need to remove default coloring of background samples
  highlight_colors = highlight_colors.map((x) => (x === "#cccccc" ? null : x));

  // Hide samples
  let data = [];
  for (let cat of dataset) data.push(cat.data);
  let plot_group_div = $("#" + target).closest(".mqc_hcplot_plotgroup");
  hide_samples(plot_group_div, samples, data);

  // Utility function to convert data for one trace into x array
  function cat_to_x(cat) {
    if (mqc_plots[target].p_active && cat["data_pct"] !== undefined) {
      // If percentage data is available, use it
      return cat["data_pct"];
    } else {
      // Otherwise use the absolute data
      return cat.data;
    }
  }

  // Render the plotly plot
  let traces = [];
  for (let cat of dataset) {
    let trace = {
      type: "bar",
      y: samples,
      x: cat_to_x(cat),
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
  // TODO: move this to multiqc_plotly.js?
  Plotly.newPlot(target, traces, layout, {
    displayModeBar: true,
    displaylogo: false,
    modeBarButtonsToRemove: [
      "lasso2d",
      "autoScale2d",
      "pan2d",
      "select2d",
      "zoom2d",
      "zoomIn2d",
      "zoomOut2d",
      "resetScale2d",
    ],
  });

  plot.replot = function () {
    // Updates plot given new underlying data
    let dataset = plot.datasets[plot.active_dataset_idx];

    let x = [];
    for (let cat of dataset) x.push(cat_to_x(cat));

    Plotly.restyle(target, "x", x);
  };
}
