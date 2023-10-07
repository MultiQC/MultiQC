// Scatter plot
function plot_scatter_plot(mqc_plots, target, dataset_idx) {
  if (mqc_plots[target] === undefined || mqc_plots[target]["plot_type"] !== "scatter") {
    return false;
  }
  var config = mqc_plots[target]["config"];
  if (dataset_idx === undefined) {
    dataset_idx = 0;
  }

  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data intact
  var data = JSON.parse(JSON.stringify(mqc_plots[target]["datasets"][dataset_idx]));

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
      if ("marker" in data[j]) {
        data[j]["marker"]["lineWidth"] = 0;
      } else {
        data[j]["marker"] = { lineWidth: 0 };
      }
      var match = false;
      $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
        if (f_text == "") {
          return true;
        }
        if (
          (window.mqc_highlight_regex_mode && data[j]["name"].match(f_text)) ||
          (!window.mqc_highlight_regex_mode && data[j]["name"].indexOf(f_text) > -1)
        ) {
          data[j]["color"] = window.mqc_highlight_f_cols[idx];
          match = true;
        }
      });
      if (!match) {
        data[j]["color"] = "rgba(100,100,100,0.2)";
      }
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

  // Make the highcharts plot
  Highcharts.chart(
    target,
    {
      chart: {
        type: "scatter",
        zoomType: "xy",
        plotBorderWidth: 1,
        height: config["square"] ? 500 : undefined,
        width: config["square"] ? 500 : undefined,
      },
      title: {
        text: config["title"],
        x: 30, // fudge to center over plot area rather than whole plot
      },
      xAxis: {
        title: {
          text: config["xlab"],
        },
        type: config["xLog"] ? "logarithmic" : "linear",
        gridLineWidth: 1,
        categories: config["categories"],
        ceiling: config["xCeiling"],
        floor: config["xFloor"],
        max: config["xmax"],
        min: config["xmin"],
        minRange: config["xMinRange"],
        allowDecimals: config["xDecimals"],
        plotBands: config["xPlotBands"],
        plotLines: config["xPlotLines"],
      },
      yAxis: {
        title: {
          text: config["ylab"],
        },
        type: config["yLog"] ? "logarithmic" : "linear",
        ceiling: config["yCeiling"],
        floor: config["yFloor"],
        max: config["ymax"],
        min: config["ymin"],
        minRange: config["yMinRange"],
        allowDecimals: config["yDecimals"],
        plotBands: config["yPlotBands"],
        plotLines: config["yPlotLines"],
      },
      plotOptions: {
        series: {
          animation: false,
          marker: {
            radius: config["marker_size"],
            lineColor: config["marker_line_colour"],
            lineWidth: config["marker_line_width"],
            states: {
              hover: {
                enabled: config["enableHover"] == undefined ? true : config["enableHover"],
                lineColor: "rgb(100,100,100)",
              },
            },
          },
          turboThreshold: config["turboThreshold"],
          enableMouseTracking: config["enableMouseTracking"],
          cursor: config["cursor"],
          point: {
            events: {
              click: config["click_func"],
            },
          },
        },
      },
      legend: {
        enabled: false,
      },
      tooltip: {
        headerFormat: "",
        pointFormat: config["pointFormat"],
        useHTML: true,
        formatter: function () {
          if (!this.point.noTooltip) {
            // Formatter function doesn't do name for some reason
            fstring = config["pointFormat"].replace("{point.name}", this.point.name);
            return Highcharts.Point.prototype.tooltipFormatter.call(this, fstring);
          }
          return false;
        },
      },
      series: [
        {
          color: config["marker_colour"],
          data: data,
        },
      ],
    },
    // Maintain aspect ratio as chart size changes
    function (this_chart) {
      if (config["square"]) {
        var resizeCh = function (chart) {
          // Extra width for legend
          var lWidth = chart.options.legend.enabled ? 30 : 0;
          // Work out new chart width, assuming needs to be narrower
          var chHeight = $(chart.renderTo).height();
          var chWidth = $(chart.renderTo).width();
          var nChHeight = chHeight;
          var nChWidth = chHeight + lWidth;
          // Chart is already too narrow, make it less tall
          if (chWidth < nChWidth) {
            nChHeight = chWidth - lWidth;
            nChWidth = chWidth;
          }
          chart.setSize(nChWidth, nChHeight);
        };
        // Resize on load
        resizeCh(this_chart);
        // Resize on graph resize
        $(this_chart.renderTo).on("mqc_plotresize", function (e) {
          resizeCh(this_chart);
        });
      }
    },
  );
}
