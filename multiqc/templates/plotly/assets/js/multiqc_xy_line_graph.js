// Basic Line Graph
function plot_xy_line_graph(plot, target, dataset_idx) {
  if (plot === undefined || plot["plot_type"] !== "xy_line") {
    return false;
  }
  var config = plot["config"];

  if (dataset_idx === undefined) {
    dataset_idx = 0;
  }
  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data intact
  var data = JSON.parse(JSON.stringify(plot["datasets"][dataset_idx]));

  // if (config["tt_label"] === undefined) {
  //   config["tt_label"] = "{point.x}: {point.y:.2f}";
  //   if (config["categories"]) {
  //     config["tt_formatter"] = function () {
  //       yval =
  //         Highcharts.numberFormat(this.y, config["tt_decimals"] == undefined ? 0 : config["tt_decimals"]) +
  //         (config["tt_suffix"] || "");
  //       return (
  //         '<div style="background-color:' +
  //         this.series.color +
  //         '; display:inline-block; height: 10px; width: 10px; border:1px solid #333;"></div> <span style="text-decoration:underline; font-weight:bold;">' +
  //         this.series.name +
  //         "</span><br><strong>" +
  //         this.key +
  //         ":</strong> " +
  //         yval
  //       );
  //     };
  //   }
  // }

  // Rename samples
  if (window.mqc_rename_f_texts.length > 0) {
    $.each(data, function (j, s) {
      $.each(window.mqc_rename_f_texts, function (idx, f_text) {
        if (window.mqc_rename_regex_mode) {
          var re = new RegExp(f_text, "g");
          data[j]["name"] = data[j]["name"].replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          data[j]["name"] = data[j]["name"].replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
    });
  }

  // Highlight samples
  if (window.mqc_highlight_f_texts.length > 0) {
    $.each(data, function (j, s) {
      $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
        if (
          (window.mqc_highlight_regex_mode && data[j]["name"].match(f_text)) ||
          (!window.mqc_highlight_regex_mode && data[j]["name"].indexOf(f_text) > -1)
        ) {
          data[j]["color"] = window.mqc_highlight_f_cols[idx];
        }
      });
    });
  }

  // Hide samples
  $("#" + target)
    .closest(".mqc_hcplot_plotgroup")
    .parent()
    .find(".samples-hidden-warning")
    .remove();
  $("#" + target)
    .closest(".mqc_hcplot_plotgroup")
    .show();
  if (window.mqc_hide_f_texts.length > 0) {
    var num_hidden = 0;
    var num_total = data.length;
    var j = data.length;
    while (j--) {
      var match = false;
      for (i = 0; i < window.mqc_hide_f_texts.length; i++) {
        var f_text = window.mqc_hide_f_texts[i];
        if (window.mqc_hide_regex_mode) {
          if (data[j]["name"].match(f_text)) {
            match = true;
          }
        } else {
          if (data[j]["name"].indexOf(f_text) > -1) {
            match = true;
          }
        }
      }
      if (window.mqc_hide_mode == "show") {
        match = !match;
      }
      if (match) {
        data.splice(j, 1);
        num_hidden += 1;
      }
    }
    // Some series hidden. Show a warning text string.
    if (num_hidden > 0) {
      var alert =
        '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> ' +
        num_hidden +
        ' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
      $("#" + target)
        .closest(".mqc_hcplot_plotgroup")
        .before(alert);
    }
    // All series hidden. Hide the graph.
    if (num_hidden == num_total) {
      $("#" + target)
        .closest(".mqc_hcplot_plotgroup")
        .hide();
      return false;
    }
  }

  // Toggle buttons for y-axis limis
  // Handler for this is at top, so doesn't get created multiple times
  if (config["ymax"] != undefined || config["ymin"] != undefined) {
    var pgroup = $("#" + target).closest(".mqc_hcplot_plotgroup");
    var wrapper = $('<div class="mqc_hcplot_yaxis_limit_toggle hidden-xs" />').prependTo(pgroup);
    wrapper.append(
      '<span class="mqc_switch_wrapper" data-ymax="' +
        config["ymax"] +
        '" data-ymin="' +
        config["ymin"] +
        '" data-target="#' +
        target +
        '">Y-Limits: <span class="mqc_switch on">on</span></span>',
    );
    wrapper.after('<div class="clearfix" />');
  }

  // Make the highcharts plot
  // Highcharts.chart(target, {
  //   chart: {
  //     type: "line",
  //     zoomType: "x",
  //   },
  //   title: {
  //     text: config["title"],
  //     x: 30, // fudge to center over plot area rather than whole plot
  //   },
  //   xAxis: {
  //     title: {
  //       text: config["xlab"],
  //     },
  //     labels: { format: config["xLabelFormat"] ? config["xLabelFormat"] : "{value}" },
  //     type: config["xLog"] ? "logarithmic" : "linear",
  //     categories: config["categories"],
  //     ceiling: config["xCeiling"],
  //     floor: config["xFloor"],
  //     max: config["xmax"],
  //     min: config["xmin"],
  //     minRange: config["xMinRange"],
  //     allowDecimals: config["xDecimals"],
  //     plotBands: config["xPlotBands"],
  //     plotLines: config["xPlotLines"],
  //   },
  //   yAxis: {
  //     title: {
  //       text: config["ylab"],
  //     },
  //     labels: { format: config["yLabelFormat"] ? config["yLabelFormat"] : "{value}" },
  //     type: config["ytype"],
  //     ceiling: config["yCeiling"],
  //     floor: config["yFloor"],
  //     max: config["ymax"],
  //     min: config["ymin"],
  //     minRange: config["yMinRange"],
  //     allowDecimals: config["yDecimals"],
  //     plotBands: config["yPlotBands"],
  //     plotLines: config["yPlotLines"],
  //   },
  //   plotOptions: {
  //     series: {
  //       marker: { enabled: false },
  //       cursor: config["cursor"],
  //       point: {
  //         events: {
  //           click: config["click_func"],
  //         },
  //       },
  //     },
  //   },
  //   legend: {
  //     enabled: false,
  //   },
  //   tooltip: {
  //     headerFormat: "",
  //     pointFormat: config["pointFormat"],
  //     formatter: config["tt_formatter"],
  //     useHTML: true,
  //   },
  //   series: data,
  // });
}
