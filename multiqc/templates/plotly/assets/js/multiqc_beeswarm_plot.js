// Beeswarm plot
function plot_beeswarm_graph(target, dataset_idx) {
  if (mqc_plots[target] === undefined || mqc_plots[target]["plot_type"] !== "beeswarm") {
    return false;
  }
  var config = mqc_plots[target]["config"];
  if (ds === undefined) {
    ds = 0;
  }

  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data in tact
  var datasets = JSON.parse(JSON.stringify(mqc_plots[target]["datasets"]));
  var samples = JSON.parse(JSON.stringify(mqc_plots[target]["samples"]));
  var categories = JSON.parse(JSON.stringify(mqc_plots[target]["categories"]));

  // Rename samples
  if (window.mqc_rename_f_texts.length > 0) {
    for (i = 0; i < samples.length; i++) {
      for (j = 0; j < samples[i].length; j++) {
        $.each(window.mqc_rename_f_texts, function (idx, f_text) {
          if (window.mqc_rename_regex_mode) {
            var re = new RegExp(f_text, "g");
            samples[i][j] = samples[i][j].replace(re, window.mqc_rename_t_texts[idx]);
          } else {
            samples[i][j] = samples[i][j].replace(f_text, window.mqc_rename_t_texts[idx]);
          }
        });
      }
    }
  }

  // Highlight samples
  var baseColour = "rgb(55,126,184)"; // Blue points by default
  var seriesColours = {};
  if (window.mqc_highlight_f_texts.length > 0) {
    baseColour = "rgb(80,80,80)"; // Grey points if no highlight
    for (i = 0; i < samples.length; i++) {
      for (j = 0; j < samples[i].length; j++) {
        $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
          if (
            (window.mqc_highlight_regex_mode && samples[i][j].match(f_text)) ||
            (!window.mqc_highlight_regex_mode && samples[i][j].indexOf(f_text) > -1)
          ) {
            seriesColours[samples[i][j]] = window.mqc_highlight_f_cols[idx];
          }
        });
      }
    }
  }

  // Hide samples
  $("#" + target)
    .closest(".hc-plot-wrapper")
    .parent()
    .find(".samples-hidden-warning")
    .remove();
  $("#" + target)
    .closest(".hc-plot-wrapper")
    .show();
  if (window.mqc_hide_f_texts.length > 0) {
    var num_hidden = 0;
    var num_total = 0;
    for (i = 0; i < samples.length; i++) {
      num_total = Math.max(num_total, samples[i].length);
      var j = samples[i].length;
      var hidden_here = 0;
      while (j--) {
        var s_name = samples[i][j];
        var match = false;
        for (k = 0; k < window.mqc_hide_f_texts.length; k++) {
          var f_text = window.mqc_hide_f_texts[k];
          if (window.mqc_hide_regex_mode) {
            if (s_name.match(f_text)) {
              match = true;
            }
          } else {
            if (s_name.indexOf(f_text) > -1) {
              match = true;
            }
          }
        }
        if (window.mqc_hide_mode == "show") {
          match = !match;
        }
        if (match) {
          samples[i].splice(j, 1);
          datasets[i].splice(j, 1);
          hidden_here += 1;
        }
      }
      num_hidden = Math.max(num_hidden, hidden_here);
    }
    // Some series hidden. Show a warning text string.
    if (num_hidden > 0) {
      var alert =
        '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> ' +
        num_hidden +
        ' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
      $("#" + target)
        .closest(".hc-plot-wrapper")
        .before(alert);
    }
    // All series hidden. Hide the graph.
    if (num_hidden == num_total) {
      $("#" + target)
        .closest(".hc-plot-wrapper")
        .hide();
      return false;
    }
  }

  // Figure out how tall to make each plot
  var ph_min = 40;
  var ph_max = 100;
  var pheight = 600 / categories.length;
  pheight = Math.min(ph_max, Math.max(ph_min, pheight));

  // Clear the loading text and add hover text placeholder
  $("#" + target).html(
    '<div class="beeswarm-hovertext"><em class="placeholder">Hover over a data point for more information</em></div><div class="beeswarm-plots"></div>',
  );
  // Resize the parent draggable div
  $("#" + target)
    .parent()
    .css("height", pheight * categories.length + 40 + "px");

  for (var i = 0; i < categories.length; i++) {
    var borderCol = categories[i]["bordercol"];
    if (borderCol == undefined) {
      borderCol = "#cccccc";
    }

    var data = datasets[i];
    var s_names = samples[i];
    if (categories[i]["namespace"] == "") {
      var label = categories[i]["title"];
      var label_long = categories[i]["description"];
    } else {
      var label = categories[i]["namespace"] + "<br/>" + categories[i]["title"];
      var label_long = categories[i]["namespace"] + ": " + categories[i]["description"];
    }
    var ttSuffix = categories[i]["suffix"];
    var decimalPlaces = categories[i]["decimalPlaces"];
    var minx = categories[i]["min"];
    var maxx = categories[i]["max"];

    // Size and spacing options
    var markerRadius = 2.5;
    var yspace = 70;
    var ysep = 10;
    if (data.length > 50) {
      markerRadius = 1.8;
      yspace = 50;
      ysep = 20;
    }
    if (data.length > 200) {
      markerRadius = 1;
      yspace = 30;
      ysep = 30;
    }

    if (maxx == undefined) {
      maxx = Math.max.apply(null, data);
    }
    if (minx == undefined) {
      minx = Math.max.apply(null, data);
    }
    var range = maxx - minx;
    var sep = range / yspace;
    // Get an array of indexes from a sorted data array
    // Leaves the data order in tact so we don't lose s_name association
    var indices = new Array(data.length);
    for (var n = 0; n < data.length; n++) {
      indices[n] = n;
    }
    indices.sort(function (a, b) {
      return data[a] < data[b] ? -1 : data[a] > data[b] ? 1 : 0;
    });
    var xydata = [];
    var last = undefined;
    var side = 1;
    for (var s_idx = 0; s_idx < indices.length; s_idx++) {
      row = indices[s_idx];
      s_name = s_names[row];
      d = data[row];
      if (Math.floor(d / sep) !== last) {
        last = Math.floor(d / sep);
        side = 1;
      } else {
        side += 1;
      }
      multiplier = side % 2 == 0 ? 1 : -1;
      var y = (Math.floor(side / 2) * multiplier) / ysep;
      // Don't let jitter get too big
      while (y > 1 || y < -1) {
        var n = Math.floor(Math.abs(y)) + 1;
        y = (Math.floor(side / 2) * multiplier) / (ysep * n);
      }
      // Get the point colour
      var thisCol = baseColour;
      if (s_name in seriesColours) {
        thisCol = seriesColours[s_name];
      }
      xydata.push({
        x: d,
        y: y,
        name: s_name,
        color: thisCol,
      });
    }

    //   $('<div class="beeswarm-plot" />')
    //     .appendTo("#" + target + " .beeswarm-plots")
    //     .css({
    //       "border-left": "2px solid " + borderCol,
    //       height: 100 / categories.length + "%",
    //     })
    //     .highcharts({
    //       chart: {
    //         type: "scatter",
    //         spacingTop: 0,
    //         marginBottom: 0,
    //         marginRight: 20,
    //         marginLeft: 180,
    //         backgroundColor: "transparent",
    //         // Horrible hacky HighCharts reflow problem.
    //         // TODO: Come back and find a better solution!
    //         events: {
    //           load: function (chart) {
    //             setTimeout(function () {
    //               chart.target.reflow();
    //             }, 200);
    //           },
    //         },
    //       },
    //       title: {
    //         text: label,
    //         align: "left",
    //         verticalAlign: "middle",
    //         y: 10,
    //         useHTML: true,
    //         style: {
    //           fontSize: "12px",
    //         },
    //       },
    //       yAxis: {
    //         title: { text: null },
    //         max: 1,
    //         min: -1,
    //         gridLineWidth: 0,
    //         title: { text: null },
    //         labels: { enabled: false },
    //         lineWidth: 0,
    //       },
    //       xAxis: {
    //         lineWidth: 0,
    //         tickWidth: 0,
    //         tickPixelInterval: 200,
    //         labels: {
    //           reserveSpace: false,
    //           y: -1 * (pheight / 2) + 5,
    //           zIndex: 1,
    //           style: {
    //             color: "#999999",
    //           },
    //         },
    //         min: minx,
    //         max: maxx,
    //       },
    //       tooltip: {
    //         valueSuffix: ttSuffix,
    //         valueDecimals: decimalPlaces,
    //         formatter: function () {
    //           var value = Highcharts.numberFormat(this.point.x, this.series.tooltipOptions.valueDecimals);
    //           var suff = this.series.tooltipOptions.valueSuffix;
    //           var ttstring =
    //             '<span style="float:right;">' +
    //             this.series.name +
    //             "</span><samp>" +
    //             this.point.name +
    //             "</samp>: &nbsp; <strong>" +
    //             value +
    //             " " +
    //             suff +
    //             "</strong>";
    //           $("#" + target + " .beeswarm-hovertext").html(ttstring);
    //           return false;
    //         },
    //       },
    //       plotOptions: {
    //         series: {
    //           name: label_long,
    //           turboThreshold: 0,
    //           marker: {
    //             radius: markerRadius,
    //             states: {
    //               hover: {
    //                 radiusPlus: 4,
    //                 lineWidthPlus: 2,
    //                 lineColor: "#333333",
    //               },
    //             },
    //           },
    //           stickyTracking: false,
    //           point: {
    //             events: {
    //               mouseOver: function (e) {
    //                 var hovName = this.name;
    //                 $("#" + target + " .beeswarm-plot").each(function () {
    //                   var plot = $(this).highcharts();
    //                   for (i = 0; i < plot.series[0].data.length; ++i) {
    //                     if (plot.series[0].data[i].name == hovName) {
    //                       plot.series[0].data[i].setState("hover");
    //                     }
    //                   }
    //                 });
    //               },
    //               mouseOut: function () {
    //                 $("#" + target + " .beeswarm-plot").each(function () {
    //                   var plot = $(this).highcharts();
    //                   for (i = 0; i < plot.series[0].data.length; ++i) {
    //                     plot.series[0].data[i].setState();
    //                   }
    //                 });
    //                 $("#" + target + " .beeswarm-hovertext").html(
    //                   '<em class="placeholder">Hover over a data point for more information</em>',
    //                 );
    //               },
    //             },
    //           },
    //         },
    //       },
    //       legend: { enabled: false },
    //       credits: { enabled: false },
    //       exporting: { enabled: false },
    //       series: [
    //         {
    //           data: xydata,
    //           // Workaround for HighCharts bug. See https://github.com/highcharts/highcharts/issues/1440
    //           marker: { states: { hover: { fillColor: {} } } },
    //         },
    //       ],
    //     });
  }
}
