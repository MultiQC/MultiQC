class HeatmapPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.xcats_samples = dump.xcats_samples;
    this.ycats_samples = dump.ycats_samples;
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    let rows = this.datasets[this.active_dataset_idx].rows;
    if (rows.length === 0) return 0; // no rows
    return rows[0].length; // no columns in a row
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let rows = this.datasets[this.active_dataset_idx].rows;
    let xcats = this.datasets[this.active_dataset_idx].xcats;
    let ycats = this.datasets[this.active_dataset_idx].ycats;

    if (this.xcats_samples) {
      let xcatsSettings = applyToolboxSettings(xcats);
      if (xcatsSettings == null) return; // All series are hidden, do not render the graph.
      // Rename and filter samples:
      xcats = xcatsSettings.filter((s) => !s.hidden).map((s) => s.name);
      rows = rows.map((row) => row.filter((val, i) => !xcatsSettings[i].hidden));
    }

    if (this.ycats_samples) {
      let ycatsSettings = applyToolboxSettings(ycats);
      if (ycatsSettings == null) return; // All series are hidden, do not render the graph.
      // Rename and filter samples:
      ycats = ycatsSettings.filter((s) => !s.hidden).map((s) => s.name);
      rows = rows.filter((row, i) => !ycatsSettings[i].hidden);
    }

    return [
      {
        z: rows,
        x: xcats,
        y: ycats,
        type: "heatmap",
      },
    ];
  }
}

// // Heatmap plot
// function plot_heatmap(plot, target, dataset_idx) {
//   if (plot === undefined || plot["plot_type"] !== "heatmap") {
//     return false;
//   }
//   var config = plot["config"];
//
//   if (config["square"] === undefined) {
//     config["square"] = true;
//   }
//   if (config["xcats_samples"] === undefined) {
//     config["xcats_samples"] = true;
//   }
//   if (config["ycats_samples"] === undefined) {
//     config["ycats_samples"] = true;
//   }
//
//   // Make a clone of the data, so that we can mess with it,
//   // while keeping the original data intact
//   var data = JSON.parse(JSON.stringify(plot["data"]));
//   var xcats = JSON.parse(JSON.stringify(plot["xcats"]));
//   var ycats = JSON.parse(JSON.stringify(plot["ycats"]));
//   // "xcats" and "ycats" are labels of columns and rows respectively
//   // data[n] has form of [x,y,value], x/y are indices of columns/rows
//
//   // Rename samples
//   if (window.mqc_rename_f_texts.length > 0) {
//     if (config["xcats_samples"]) {
//       for (i = 0; i < xcats.length; i++) {
//         $.each(window.mqc_rename_f_texts, function (idx, f_text) {
//           if (window.mqc_rename_regex_mode) {
//             var re = new RegExp(f_text, "g");
//             xcats[i] = xcats[i].replace(re, window.mqc_rename_t_texts[idx]);
//           } else {
//             xcats[i] = xcats[i].replace(f_text, window.mqc_rename_t_texts[idx]);
//           }
//         });
//       }
//     }
//     if (config["ycats_samples"]) {
//       for (i = 0; i < ycats.length; i++) {
//         $.each(window.mqc_rename_f_texts, function (idx, f_text) {
//           if (window.mqc_rename_regex_mode) {
//             var re = new RegExp(f_text, "g");
//             ycats[i] = ycats[i].replace(re, window.mqc_rename_t_texts[idx]);
//           } else {
//             ycats[i] = ycats[i].replace(f_text, window.mqc_rename_t_texts[idx]);
//           }
//         });
//       }
//     }
//   }
//
//   // Sort samples by highlight
//   $(".mqc_heatmap_sortHighlight").attr("disabled", false);
//   if (config["sortHighlights"] == true) {
//     if (window.mqc_highlight_f_texts.length > 0) {
//       // Collect the highlighting indices
//       var xcat_hl = Array();
//       var ycat_hl = Array();
//       for (i = 0; i < xcats.length; i++) {
//         $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
//           if (f_text == "") {
//             xcat_hl[i] = 0;
//           } else if (
//             (window.mqc_highlight_regex_mode && xcats[i].match(f_text)) ||
//             (!window.mqc_highlight_regex_mode && xcats[i].indexOf(f_text) > -1)
//           ) {
//             xcat_hl[i] = window.mqc_highlight_f_texts.length - idx;
//           }
//         });
//       }
//       for (i = 0; i < ycats.length; i++) {
//         $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
//           if (f_text == "") {
//             ycat_hl[i] = 0;
//           } else if (
//             (window.mqc_highlight_regex_mode && ycats[i].match(f_text)) ||
//             (!window.mqc_highlight_regex_mode && ycats[i].indexOf(f_text) > -1)
//           ) {
//             ycat_hl[i] = window.mqc_highlight_f_texts.length - idx;
//           }
//         });
//       }
//       // Reshape the data - needs deepcopy as indexes are updated
//       var newdata = JSON.parse(JSON.stringify(plot["data"]));
//       var new_xcats = [],
//         new_ycats = [];
//       var xidx = 0,
//         yidx = 0;
//       for (hl = window.mqc_highlight_f_texts.length; hl >= 0; hl--) {
//         if (config["xcats_samples"]) {
//           for (i = 0; i < xcats.length; i++) {
//             if (xcat_hl[i] == hl) {
//               new_xcats.push(xcats[i]);
//               for (j = 0; j < data.length; j++) {
//                 // data[j] element is [x,y,val], get "x"
//                 if (data[j][0] == i) {
//                   newdata[j][0] = xidx;
//                 }
//               }
//               xidx += 1;
//             }
//           }
//         }
//         if (config["ycats_samples"]) {
//           for (i = 0; i < ycats.length; i++) {
//             if (ycat_hl[i] == hl) {
//               new_ycats.push(ycats[i]);
//               for (j = 0; j < data.length; j++) {
//                 // data[j] element is [x,y,val], get "y"
//                 if (data[j][1] == i) {
//                   newdata[j][1] = yidx;
//                 }
//               }
//               yidx += 1;
//             }
//           }
//         }
//       }
//       data = newdata;
//       if (config["xcats_samples"]) {
//         xcats = new_xcats;
//       }
//       if (config["ycats_samples"]) {
//         ycats = new_ycats;
//       }
//     }
//   }
//
//   // Hide samples
//   var num_total = Math.max(xcats.length, ycats.length);
//   $("#" + target)
//     .closest(".hc-plot-wrapper")
//     .parent()
//     .find(".samples-hidden-warning")
//     .remove();
//   $("#" + target)
//     .closest(".hc-plot-wrapper")
//     .show();
//   if (window.mqc_hide_f_texts.length > 0) {
//     var remove = Array();
//     if (config["xcats_samples"]) {
//       var i = xcats.length;
//       var xhidden = 0;
//       // iterate over x-categories (columns)
//       while (i--) {
//         var match = false;
//         for (j = 0; j < window.mqc_hide_f_texts.length; j++) {
//           var f_text = window.mqc_hide_f_texts[j];
//           if (window.mqc_hide_regex_mode) {
//             if (xcats[i].match(f_text)) {
//               match = true;
//             }
//           } else {
//             if (xcats[i].indexOf(f_text) > -1) {
//               match = true;
//             }
//           }
//         }
//         if (window.mqc_hide_mode == "show") {
//           match = !match;
//         }
//         // modify data if "i" is match for "hiding":
//         // mark elements from "i"-th column (x) for removal,
//         // shift "x" of elements from "i+" columns to the left
//         if (match) {
//           xcats.splice(i, 1);
//           for (n = 0; n < data.length; n++) {
//             // data[n] element is [x,y,val], get "x"
//             let x = data[n][0];
//             if (x == i) {
//               remove.push(n);
//             } else if (x > i) {
//               data[n][0] = x - 1;
//             }
//           }
//           xhidden += 1;
//         }
//       }
//     }
//     if (config["ycats_samples"]) {
//       var i = ycats.length;
//       var yhidden = 0;
//       // iterate over y-categories (rows)
//       while (i--) {
//         var match = false;
//         for (j = 0; j < window.mqc_hide_f_texts.length; j++) {
//           var f_text = window.mqc_hide_f_texts[j];
//           if (window.mqc_hide_regex_mode) {
//             if (ycats[i].match(f_text)) {
//               match = true;
//             }
//           } else {
//             if (ycats[i].indexOf(f_text) > -1) {
//               match = true;
//             }
//           }
//         }
//         if (window.mqc_hide_mode == "show") {
//           match = !match;
//         }
//         // modify data if "i" is match for "hiding"
//         // mark elements from "i"-th row (y) for removal,
//         // shift "y" of elements from "i+" rows up
//         if (match) {
//           ycats.splice(i, 1);
//           for (n = 0; n < data.length; n++) {
//             // data[n] element is [x,y,val], get "y"
//             let y = data[n][1];
//             if (y == i) {
//               if (remove.indexOf(n) < 0) {
//                 remove.push(n);
//               }
//             } else if (y > i) {
//               data[n][1] = y - 1;
//             }
//           }
//           yhidden += 1;
//         }
//       }
//     }
//     // Remove the data values that matched
//     remove = remove.sort(function (a, b) {
//       return a - b;
//     }); // Sorts alphabetically by default, even with integers
//     var r = remove.length;
//     while (r--) {
//       data.splice(remove[r], 1);
//     }
//     // Report / hide the plot if we're hiding stuff
//     var num_hidden = Math.max(xhidden, yhidden);
//     // Some series hidden. Show a warning text string.
//     if (num_hidden > 0) {
//       var alert =
//         '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> ' +
//         num_hidden +
//         ' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
//       $("#" + target)
//         .closest(".hc-plot-wrapper")
//         .before(alert);
//     }
//     // All series hidden. Hide the graph.
//     if (num_hidden >= num_total) {
//       $("#" + target)
//         .closest(".hc-plot-wrapper")
//         .hide();
//       return false;
//     }
//   }
//
//   // Highlight samples - do this last as we convert numerical arrays to associative
//   if (window.mqc_highlight_f_texts.length > 0) {
//     $(".mqc_heatmap_sortHighlight").attr("disabled", false);
//     var highlight_cells = Array();
//     if (config["ycats_samples"]) {
//       for (i = 0; i < xcats.length; i++) {
//         $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
//           if (f_text == "") {
//             return true;
//           }
//           if (
//             (window.mqc_highlight_regex_mode && xcats[i].match(f_text)) ||
//             (!window.mqc_highlight_regex_mode && xcats[i].indexOf(f_text) > -1)
//           ) {
//             for (n = 0; n < data.length; n++) {
//               highlight_cells[idx] =
//                 typeof highlight_cells[idx] != "undefined" && highlight_cells[idx] instanceof Array
//                   ? highlight_cells[idx]
//                   : [];
//               // data[n] element is [x,y,val], get "x"
//               if (data[n][0] == i) {
//                 highlight_cells[idx].push(n);
//               }
//             }
//           }
//         });
//       }
//     }
//     if (config["ycats_samples"]) {
//       for (i = 0; i < ycats.length; i++) {
//         $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
//           if (f_text == "") {
//             return true;
//           }
//           if (
//             (window.mqc_highlight_regex_mode && ycats[i].match(f_text)) ||
//             (!window.mqc_highlight_regex_mode && ycats[i].indexOf(f_text) > -1)
//           ) {
//             for (n = 0; n < data.length; n++) {
//               // data[n] element is [x,y,val], get "y"
//               if (data[n][1] == i) {
//                 highlight_cells[idx] =
//                   typeof highlight_cells[idx] != "undefined" && highlight_cells[idx] instanceof Array
//                     ? highlight_cells[idx]
//                     : [];
//                 if (highlight_cells[idx].indexOf(n) < 0) {
//                   highlight_cells[idx].push(n);
//                 }
//               }
//             }
//           }
//         });
//       }
//     }
//     // Give highlighted cells a border
//     for (var idx in highlight_cells) {
//       var hl = highlight_cells[idx];
//       hl = hl.sort(function (a, b) {
//         return a - b;
//       }); // Sorts alphabetically by default, even with integers
//       var h = hl.length;
//       while (h--) {
//         var i = hl[h];
//         data[i] = {
//           // data[i] element is [x,y,val]
//           x: data[i][0] === undefined ? data[i]["x"] : data[i][0],
//           y: data[i][1] === undefined ? data[i]["y"] : data[i][1],
//           value: data[i][2] === undefined ? data[i]["value"] : data[i][2],
//           borderWidth: 2,
//           borderColor: window.mqc_highlight_f_cols[idx],
//         };
//       }
//     }
//   } else {
//     $(".mqc_heatmap_sortHighlight").attr("disabled", true);
//   }
//
//   // Build scale
//   if (config["colstops"] === undefined) {
//     config["colstops"] = [
//       [0, "#313695"],
//       [0.1, "#4575b4"],
//       [0.2, "#74add1"],
//       [0.3, "#abd9e9"],
//       [0.4, "#e0f3f8"],
//       [0.5, "#ffffbf"],
//       [0.6, "#fee090"],
//       [0.7, "#fdae61"],
//       [0.8, "#f46d43"],
//       [0.9, "#d73027"],
//       [1, "#a50026"],
//     ];
//   }
//   if (config["reverseColors"] === undefined) {
//     config["reverseColors"] = false;
//   }
//   if (config["decimalPlaces"] === undefined) {
//     config["decimalPlaces"] = 2;
//   }
//   if (config["legend"] === undefined) {
//     config["legend"] = true;
//   }
//   if (config["borderWidth"] === undefined) {
//     config["borderWidth"] = 0;
//   }
//   var datalabels = config["datalabels"];
//   if (datalabels === undefined) {
//     if (data.length < 20) {
//       datalabels = true;
//     } else {
//       datalabels = false;
//     }
//   }
//   // Clone the colstops before we mess around with them
//   var colstops = JSON.parse(JSON.stringify(config["colstops"]));
//   // Reverse the colour scale if the axis is reversed
//   if (config["reverseColors"]) {
//     for (var i = 0; i < colstops.length; i++) {
//       colstops[i][0] = 1 - colstops[i][0];
//     }
//     colstops.reverse();
//   }
//
//   // // Make the highcharts plot
//   // Highcharts.chart(
//   //   target,
//   //   {
//   //     chart: {
//   //       type: "heatmap",
//   //       zoomType: "xy",
//   //       height: config["square"] ? 500 : undefined,
//   //       width: config["square"] ? 530 : undefined,
//   //       marginTop: config["title"] ? 60 : 50,
//   //     },
//   //     plotOptions: {
//   //       series: {
//   //         states: {
//   //           hover: {
//   //             borderWidth: 2,
//   //             borderColor: "red",
//   //           },
//   //         },
//   //       },
//   //     },
//   //     title: {
//   //       text: config["title"],
//   //     },
//   //     xAxis: {
//   //       endOnTick: false,
//   //       maxPadding: 0,
//   //       categories: xcats,
//   //       title: { enabled: true, text: config["xTitle"] },
//   //       labels: {
//   //         formatter: function () {
//   //           try {
//   //             return this.value.substr(0, 20);
//   //           } catch (err) {
//   //             return this.value;
//   //           }
//   //         },
//   //       },
//   //     },
//   //     yAxis: {
//   //       endOnTick: false,
//   //       maxPadding: 0,
//   //       categories: ycats,
//   //       reversed: true,
//   //       opposite: true,
//   //       title: config["yTitle"],
//   //       labels: {
//   //         formatter: function () {
//   //           try {
//   //             return this.value.substr(0, 20);
//   //           } catch (err) {
//   //             return this.value;
//   //           }
//   //         },
//   //       },
//   //     },
//   //     colorAxis: {
//   //       reversed: config["reverseColors"],
//   //       stops: colstops,
//   //       min: config["min"],
//   //       max: config["max"],
//   //     },
//   //     legend: {
//   //       align: "right",
//   //       layout: "vertical",
//   //       margin: 0,
//   //       verticalAlign: "top",
//   //       y: 25,
//   //       symbolHeight: 280,
//   //       enabled: config["legend"],
//   //     },
//   //     tooltip: {
//   //       useHTML: true,
//   //       formatter: function () {
//   //         return (
//   //           'X: <span style="font-weight:bold; font-family:monospace;">' +
//   //           this.series.xAxis.categories[this.point.x] +
//   //           "</span><br>" +
//   //           'Y: <span style="font-weight:bold; font-family:monospace;">' +
//   //           this.series.yAxis.categories[this.point.y] +
//   //           "</span><br>" +
//   //           '<div style="background-color:' +
//   //           this.point.color +
//   //           '; display:inline-block; height: 10px; width: 10px; border:1px solid #333;"></div> ' +
//   //           '<span style="font-weight: bold; text-decoration:underline;">' +
//   //           Highcharts.numberFormat(this.point.value, config["decimalPlaces"]) +
//   //           "</span>"
//   //         );
//   //       },
//   //     },
//   //     series: [
//   //       {
//   //         turboThreshold: 0,
//   //         borderWidth: config["borderWidth"],
//   //         data: data,
//   //         dataLabels: {
//   //           enabled: datalabels,
//   //           format: "{point.value:." + config["decimalPlaces"] + "f}",
//   //           color: config["datalabel_colour"],
//   //         },
//   //       },
//   //     ],
//   //   },
//   //   function (this_chart) {
//   //     // Maintain aspect ratio as chart size changes
//   //     if (config["square"]) {
//   //       var resizeCh = function (chart) {
//   //         // Extra width for legend
//   //         var lWidth = chart.options.legend.enabled ? 30 : 0;
//   //         // Work out new chart width, assuming needs to be narrower
//   //         var chHeight = $(chart.renderTo).height();
//   //         var chWidth = $(chart.renderTo).width();
//   //         var nChHeight = chHeight;
//   //         var nChWidth = chHeight + lWidth;
//   //         // Chart is already too narrow, make it less tall
//   //         if (chWidth < nChWidth) {
//   //           nChHeight = chWidth - lWidth;
//   //           nChWidth = chWidth;
//   //         }
//   //         chart.setSize(nChWidth, nChHeight);
//   //       };
//   //       // Resize on load
//   //       resizeCh(this_chart);
//   //       // Resize on graph resize
//   //       $(this_chart.renderTo).on("mqc_plotresize", function (e) {
//   //         try {
//   //           resizeCh(this_chart);
//   //         } catch (e) {
//   //           plot_heatmap($(this).attr("id"));
//   //         }
//   //       });
//   //     }
//   //   },
//   // );
//
//   // Listeners for range slider
//   $(".mqc_hcplot_range_sliders input").on("keyup change input", function () {
//     target = $(this).data("target");
//     minmax = $(this).data("minmax");
//     var chart = $("#" + target).highcharts();
//     if (minmax == "min") {
//       chart.colorAxis[0].update({ min: $(this).val() });
//     }
//     if (minmax == "max") {
//       chart.colorAxis[0].update({ max: $(this).val() });
//     }
//     $("#" + target + "_range_slider_" + minmax + ", #" + target + "_range_slider_" + minmax + "_txt").val(
//       $(this).val(),
//     );
//   });
// }
