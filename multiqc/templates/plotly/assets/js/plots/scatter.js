class ScatterPlot extends Plot {
  constructor(data) {
    super(data);
    this.categories = data.categories;
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.active_dataset_idx].points; // no data points in a dataset
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let points = this.datasets[this.active_dataset_idx].points;
    if (points.length === 0) return [];

    let samples = points.map((point) => point.name);
    let sampleSettings = applyToolboxSettings(samples);
    if (sampleSettings == null) return; // All series are hidden, do not render the graph
    points = points.map((point, idx) => {
      point.name = sampleSettings[idx].name ?? point.name;
      point.highlight = sampleSettings[idx].highlight;
      if (!sampleSettings[idx].hidden) return point;
    });

    // Reorder points so highlighted points are on top
    let highlighted = points.filter((p) => p.highlight);
    let nonHighlighted = points.filter((p) => !p.highlight);
    points = nonHighlighted.concat(highlighted);

    return points.map((point) => {
      let x = point.x;
      if (this.categories && Number.isInteger(x) && x < this.categories.length) x = this.categories[x];

      let params = JSON.parse(JSON.stringify(this.trace_params)); // deep copy
      params.marker.size = point["marker_size"] ?? params.marker.size;
      params.marker.line = {
        width: point["marker_line_width"] ?? params.marker.line.width,
      };
      params.marker.opacity = point["opacity"] ?? params.marker.opacity;
      params.marker.color = point["color"] ?? params.marker.color;
      if (highlighted.length > 0) params.marker.color = point.highlight ?? "#cccccc";

      return {
        type: "scatter",
        x: [x],
        y: [point.y],
        name: point.name,
        text: [point.annotation ?? point.name],
        ...params,
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
