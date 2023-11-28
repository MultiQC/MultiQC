class LinePlot extends Plot {
  constructor(target, data) {
    super(target, data);
    this.datasets = data.datasets;
  }

  active_dataset_size() {
    if (this.datasets.length === 0) return 0;
    return this.datasets[this.active_dataset_idx].length;
  }

  // Constructs and returns traces for the Plotly plot
  build_traces() {
    let dataset = this.datasets[this.active_dataset_idx];

    // Samples for active dataset
    let samples = [];
    for (let sdata of dataset) samples.push(sdata.name);

    // Rename samples
    rename_samples(samples, function (sample_idx, new_name) {
      dataset[sample_idx]["name"] = new_name;
    });

    // Hide samples
    let plot_group_div = $("#" + this.target).closest(".mqc_hcplot_plotgroup");
    let idx_to_hide = hide_samples(plot_group_div, samples);
    if (idx_to_hide.length === samples.length)
      // All series hidden. Hide the graph.
      return;
    let visible_samples = samples.filter((x, i) => !idx_to_hide.includes(i));
    let visible_dataset = dataset.filter((x, i) => !idx_to_hide.includes(i));

    // Highlight samples
    let highlight_colors = get_highlight_colors(visible_samples);

    // Toggle buttons for Y-axis limis
    // Handler for this is at top, so doesn't get created multiple times
    let ymax_is_set = this.pconfig.ymax !== "undefined" && this.pconfig.ymax !== null;
    let ymin_is_set = this.pconfig.ymin !== "undefined" && this.pconfig.ymin !== null;
    // Only create if there is a y-axis limit
    if (ymax_is_set || (ymin_is_set && this.pconfig.ymin !== 0)) {
      let wrapper = $('<div class="mqc_hcplot_yaxis_limit_toggle hidden-xs" />').prependTo(plot_group_div);
      wrapper.append(
        '<span class="mqc_switch_wrapper"' +
          '"' +
          ' data-ymax="' +
          this.pconfig.ymax +
          '"' +
          ' data-ymin="' +
          this.pconfig.ymin +
          '"' +
          ' data-target="' +
          this.target +
          '">Y-Limits: <span class="mqc_switch on">on</span></span>',
      );
      wrapper.after('<div class="clearfix" />');
    }

    return visible_dataset.map((sdata, si) => {
      // let [x, y] = this.sdata_to_xy(dataset[i]);
      let x, y;
      if (sdata.data.length > 0 && Array.isArray(sdata.data[0])) {
        x = sdata.data.map((x) => x[0]);
        y = sdata.data.map((x) => x[1]);
      } else {
        x = [...Array(sdata.data.length).keys()];
        y = sdata.data;
      }
      return {
        type: "scatter",
        x: x,
        y: y,
        name: sdata.name,
        orientation: "h",
        mode: "lines",
        line: {
          color: highlight_colors[si] || sdata.color,
        },
        // hovertemplate: pconfig["tt_label"]
      };
    });
  }
}
//   // Utility function to convert data for one trace into x and y arrays
//   sdata_to_xy(sdata) {
//     let x, y;
//     if (sdata.data.length > 0 && Array.isArray(sdata.data[0])) {
//       x = sdata.data.map((x) => x[0]);
//       y = sdata.data.map((x) => x[1]);
//     } else {
//       x = [...Array(sdata.data.length).keys()];
//       y = sdata.data;
//     }
//     return [x, y];
//   }
// }

//   // Function to re-render plot with new data  // TODO: easier to redraw everything?
//   replot() {
//     // Updates plot given new underlying data
//     let dataset = this.datasets[this.active_dataset_idx];
//
//     let xs = [];
//     let ys = [];
//     for (let sdata of dataset) {
//       let [x, y] = this.sdata_to_xy(sdata);
//       xs.push(x);
//       ys.push(y);
//     }
//     Plotly.restyle(this.target, "x", xs);
//     Plotly.restyle(this.target, "y", ys);
//     // For some reason range breaks when restyling with new x and y
//     // Figure out why?
//     // let ymax_is_set = this.pconfig.ymax !== "undefined" && this.pconfig.ymax !== null;
//     // let ymin_is_set = this.pconfig.ymin !== "undefined" && this.pconfig.ymin !== null;
//     // let xmax_is_set = this.pconfig.xmax !== "undefined" && this.pconfig.xmax !== null;
//     // let xmin_is_set = this.pconfig.xmin !== "undefined" && this.pconfig.xmin !== null;
//     // if (ymin_is_set) Plotly.relayout(target, "yaxis.range[0]", this.pconfig.ymin);
//     // if (ymax_is_set) Plotly.relayout(target, "yaxis.range[0]", this.pconfig.ymax);
//     // if (xmin_is_set) Plotly.relayout(target, "xaxis.range[0]", this.pconfig.xmin);
//     // if (xmax_is_set) Plotly.relayout(target, "xaxis.range[0]", this.pconfig.xmax);
//   }
// }
