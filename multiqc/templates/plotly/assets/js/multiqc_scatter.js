class ScatterPlot extends Plot {
  constructor(data) {
    super(data);
    this.categories = data.categories;
    this.default_marker = data.default_marker;
  }

  activeDatasetSamples() {
    if (this.datasets.length === 0) return [];
    let dataset = this.datasets[this.active_dataset_idx];
    return dataset.map((element) => element.name);
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let dataset = this.datasets[this.active_dataset_idx];

    let samples = applyToolboxSettings(this.activeDatasetSamples());

    // All series are hidden, do not render the graph.
    if (samples == null) return;

    // Rename samples
    dataset.map((element, si) => (element.name = samples[si].name));

    // Filter out hidden samples
    let visibleSamples = samples.filter((s) => !s.hidden);
    let visibleDataset = dataset.filter((element, si) => !samples[si].hidden);

    return visibleDataset.map((element, si) => {
      let x = element.x;
      if (this.categories && Number.isInteger(x) && x < this.categories.length) x = this.categories[x];

      // Shallow copy of default marker
      let marker = Object.assign({}, this.default_marker);
      marker.size = element["marker_size"] ?? marker.size;
      marker.line.width = element["marker_line_width"] ?? marker.line.width;
      marker.color = visibleSamples[si].highlight ?? element["color"] ?? marker.color;
      marker.opacity = element["opacity"] ?? marker.opacity;

      return {
        type: "scatter",
        x: [x],
        y: [element.y],
        name: element.name,
        mode: "markers",
        marker: marker,
        // marker: {
        //   // TODO: default size, width, color, etc are repetead here and in flat plotting python code
        //   // we need to put these values into a plot config and reuse them here
        //   // we can also sync python class and JS class by building the same pconfig and putting
        //   // those values there, insted of saving __dict__ in python and then parsing it in JS as pconfig
        //   size: element["marker_size"] ?? this.layout.marker.size,
        //   line: {
        //     width: element["marker_line_width"] ?? marker["line"].width,
        //   },
        //   color: samples[element.name].highlight ?? element["color"] ?? marker["color"],
        //   opacity: element["opacity"] ?? marker["opacity"],
        // },
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
