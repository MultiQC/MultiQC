// Scatter plot
class ScatterPlot extends Plot {
  constructor(target, data) {
    super(target, data);
  }

  activeDatasetSize() {}

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let dataset = this.datasets[this.active_dataset_idx];

    // Samples for active dataset
    let samples = [];
    for (let sdata of dataset) samples.push(sdata.name);

    // Rename samples
    rename_samples(dataset.samples, function (sample_idx, new_name) {
      dataset.samples[sample_idx] = new_name;
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

    // // Highlight samples
    // if (window.mqc_highlight_f_texts.length > 0) {
    //   $.each(data, function (j, s) {
    //     if ("marker" in data[j]) {
    //       data[j]["marker"]["lineWidth"] = 0;
    //     } else {
    //       data[j]["marker"] = { lineWidth: 0 };
    //     }
    //     var match = false;
    //     $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
    //       if (f_text == "") {
    //         return true;
    //       }
    //       if (
    //         (window.mqc_highlight_regex_mode && data[j]["name"].match(f_text)) ||
    //         (!window.mqc_highlight_regex_mode && data[j]["name"].indexOf(f_text) > -1)
    //       ) {
    //         data[j]["color"] = window.mqc_highlight_f_cols[idx];
    //         match = true;
    //       }
    //     });
    //     if (!match) {
    //       data[j]["color"] = "rgba(100,100,100,0.2)";
    //     }
    //   });
    // }

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
        mode: "markers",
        marker: {
          color: highlight_colors[si] || sdata.color,
        },
        // hovertemplate: pconfig["tt_label"]
      };
    });
  }
}
// // Make the highcharts plot
// Highcharts.chart(
//   target,
//   {
//     chart: {
//       type: "scatter",
//       zoomType: "xy",
//       plotBorderWidth: 1,
//       height: config["square"] ? 500 : undefined,
//       width: config["square"] ? 500 : undefined,
//     },
//     title: {
//       text: config["title"],
//       x: 30, // fudge to center over plot area rather than whole plot
//     },
//     xAxis: {
//       title: {
//         text: config["xlab"],
//       },
//       type: config["xLog"] ? "logarithmic" : "linear",
//       gridLineWidth: 1,
//       categories: config["categories"],
//       ceiling: config["xCeiling"],
//       floor: config["xFloor"],
//       max: config["xmax"],
//       min: config["xmin"],
//       minRange: config["xMinRange"],
//       allowDecimals: config["xDecimals"],
//       plotBands: config["xPlotBands"],
//       plotLines: config["xPlotLines"],
//     },
//     yAxis: {
//       title: {
//         text: config["ylab"],
//       },
//       type: config["yLog"] ? "logarithmic" : "linear",
//       ceiling: config["yCeiling"],
//       floor: config["yFloor"],
//       max: config["ymax"],
//       min: config["ymin"],
//       minRange: config["yMinRange"],
//       allowDecimals: config["yDecimals"],
//       plotBands: config["yPlotBands"],
//       plotLines: config["yPlotLines"],
//     },
//     plotOptions: {
//       series: {
//         animation: false,
//         marker: {
//           radius: config["marker_size"],
//           lineColor: config["marker_line_colour"],
//           lineWidth: config["marker_line_width"],
//           states: {
//             hover: {
//               enabled: config["enableHover"] == undefined ? true : config["enableHover"],
//               lineColor: "rgb(100,100,100)",
//             },
//           },
//         },
//         turboThreshold: config["turboThreshold"],
//         enableMouseTracking: config["enableMouseTracking"],
//         cursor: config["cursor"],
//         point: {
//           events: {
//             click: config["click_func"],
//           },
//         },
//       },
//     },
//     legend: {
//       enabled: false,
//     },
//     tooltip: {
//       headerFormat: "",
//       pointFormat: config["pointFormat"],
//       useHTML: true,
//       formatter: function () {
//         if (!this.point.noTooltip) {
//           // Formatter function doesn't do name for some reason
//           fstring = config["pointFormat"].replace("{point.name}", this.point.name);
//           return Highcharts.Point.prototype.tooltipFormatter.call(this, fstring);
//         }
//         return false;
//       },
//     },
//     series: [
//       {
//         color: config["marker_colour"],
//         data: data,
//       },
//     ],
//   },
//   // Maintain aspect ratio as chart size changes
//   function (this_chart) {
//     if (config["square"]) {
//       var resizeCh = function (chart) {
//         // Extra width for legend
//         var lWidth = chart.options.legend.enabled ? 30 : 0;
//         // Work out new chart width, assuming needs to be narrower
//         var chHeight = $(chart.renderTo).height();
//         var chWidth = $(chart.renderTo).width();
//         var nChHeight = chHeight;
//         var nChWidth = chHeight + lWidth;
//         // Chart is already too narrow, make it less tall
//         if (chWidth < nChWidth) {
//           nChHeight = chWidth - lWidth;
//           nChWidth = chWidth;
//         }
//         chart.setSize(nChWidth, nChHeight);
//       };
//       // Resize on load
//       resizeCh(this_chart);
//       // Resize on graph resize
//       $(this_chart.renderTo).on("mqc_plotresize", function (e) {
//         resizeCh(this_chart);
//       });
//     }
//   },
// );
