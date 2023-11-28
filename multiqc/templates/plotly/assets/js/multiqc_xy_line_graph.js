// Basic Line Graph
function plot_xy_line_graph(plot, target, dataset_idx) {
  if (plot === undefined || plot["plot_type"] !== "xy_line") return false;

  if (dataset_idx === undefined) dataset_idx = 0;

  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data intact
  let data = JSON.parse(JSON.stringify(plot["datasets"][dataset_idx]));
  let layout = JSON.parse(JSON.stringify(plot["layout"]));
  let pconfig = JSON.parse(JSON.stringify(plot["pconfig"]));

  let samples = [];
  for (let sdata of data) samples.push(sdata.name);

  // Rename samples
  rename_samples(samples, function (sample_idx, new_name) {
    data[sample_idx]["name"] = new_name;
  });

  // Highlight samples
  let highlight_colors = get_highlight_colors(samples);

  // Hide samples
  let plot_group_div = $("#" + target).closest(".mqc_hcplot_plotgroup");
  hide_samples(plot_group_div, samples, data);

  // Toggle buttons for Y-axis limis
  // Handler for this is at top, so doesn't get created multiple times
  let ymax_is_set = pconfig.ymax !== "undefined" && pconfig.ymax !== null;
  let ymin_is_set = pconfig.ymin !== "undefined" && pconfig.ymin !== null;
  // Only create if there is a y-axis limit
  if (ymax_is_set || (ymin_is_set && pconfig.ymin !== 0)) {
    let wrapper = $('<div class="mqc_hcplot_yaxis_limit_toggle hidden-xs" />').prependTo(plot_group_div);
    wrapper.append(
      '<span class="mqc_switch_wrapper"' +
        '"' +
        ' data-ymax="' +
        pconfig.ymax +
        '"' +
        ' data-ymin="' +
        pconfig.ymin +
        '"' +
        ' data-target="' +
        target +
        '">Y-Limits: <span class="mqc_switch on">on</span></span>',
    );
    wrapper.after('<div class="clearfix" />');
  }

  function sdata_to_xy(sdata) {
    // Utility function to convert data for one trace into x and y arrays
    let x, y;
    if (sdata.data.length > 0 && Array.isArray(sdata.data[0])) {
      x = sdata.data.map((x) => x[0]);
      y = sdata.data.map((x) => x[1]);
    } else {
      x = [...Array(sdata.data.length).keys()];
      y = sdata.data;
    }
    return [x, y];
  }

  let traces = [];
  for (let i = 0; i < data.length; i++) {
    let [x, y] = sdata_to_xy(data[i]);
    let trace = {
      type: "scatter",
      x: x,
      y: y,
      name: data[i].name,
      orientation: "h",
      mode: "lines",
      line: {
        color: highlight_colors[i] || data[i].color,
      },
      // hovertemplate: pconfig["tt_label"]
    };
    traces.push(trace);
  }
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

  // Function to re-render plot with new data  // TODO: easier to redraw everything?
  plot.replot = function () {
    // Updates plot given new underlying data
    let dataset = plot.datasets[plot.active_dataset_idx];

    let xs = [];
    let ys = [];
    for (let sdata of dataset) {
      let [x, y] = sdata_to_xy(sdata);
      xs.push(x);
      ys.push(y);
    }
    Plotly.restyle(target, "x", xs);
    Plotly.restyle(target, "y", ys);
    // For some reason range breaks when restyling with new x and y
    // Figure out why?
    let ymax_is_set = pconfig.ymax !== "undefined" && pconfig.ymax !== null;
    let ymin_is_set = pconfig.ymin !== "undefined" && pconfig.ymin !== null;
    let xmax_is_set = pconfig.xmax !== "undefined" && pconfig.xmax !== null;
    let xmin_is_set = pconfig.xmin !== "undefined" && pconfig.xmin !== null;
    if (ymin_is_set) Plotly.relayout(target, "yaxis.range[0]", pconfig.ymin);
    if (ymax_is_set) Plotly.relayout(target, "yaxis.range[0]", pconfig.ymax);
    if (xmin_is_set) Plotly.relayout(target, "xaxis.range[0]", pconfig.xmin);
    if (xmax_is_set) Plotly.relayout(target, "xaxis.range[0]", pconfig.xmax);
  };
}
