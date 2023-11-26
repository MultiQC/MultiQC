// Basic Line Graph
function plot_xy_line_graph(plot, target, dataset_idx) {
  if (plot === undefined || plot["plot_type"] !== "xy_line") {
    return false;
  }
  if (dataset_idx === undefined) {
    dataset_idx = 0;
  }
  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data intact
  let data = JSON.parse(JSON.stringify(plot["datasets"][dataset_idx]));
  let layout = JSON.parse(JSON.stringify(plot["layout"]));
  let settings = JSON.parse(JSON.stringify(plot["settings"]));

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
    $.each(data, function (sample_idx, s_name) {
      $.each(window.mqc_rename_f_texts, function (idx, f_text) {
        if (window.mqc_rename_regex_mode) {
          const re = new RegExp(f_text, "g");
          data[sample_idx]["name"] = data[sample_idx]["name"].replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          data[sample_idx]["name"] = data[sample_idx]["name"].replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
    });
  }

  // Highlight samples
  let highlight_colors = [];
  if (window.mqc_highlight_f_texts.length > 0) {
    $.each(data, function (sample_idx, s_name) {
      highlight_colors[sample_idx] = null;
      $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
        if (f_text === "") {
          return true;
        } // skip blanks
        if (
          (window.mqc_highlight_regex_mode && s_name.match(f_text)) ||
          (!window.mqc_highlight_regex_mode && s_name.indexOf(f_text) > -1)
        ) {
          // Make the data point in each series with this index have a border colour
          highlight_colors[sample_idx] = window.mqc_highlight_f_cols[idx];
        }
      });
    });
  }

  // Hide samples
  let plot_group_div = $("#" + target).closest(".mqc_hcplot_plotgroup");
  plot_group_div.parent().find(".samples-hidden-warning").remove();
  plot_group_div.show();
  if (window.mqc_hide_f_texts.length > 0) {
    let num_hidden = 0;
    let num_total = data.length;
    let j = data.length;
    while (j--) {
      let match = false;
      for (let i = 0; i < window.mqc_hide_f_texts.length; i++) {
        const f_text = window.mqc_hide_f_texts[i];
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
      if (window.mqc_hide_mode === "show") {
        match = !match;
      }
      if (match) {
        data.splice(j, 1);
        num_hidden += 1;
      }
    }
    // Some series hidden. Show a warning text string.
    if (num_hidden > 0) {
      const alert =
        '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> ' +
        num_hidden +
        ' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
      plot_group_div.before(alert);
    }
    // All series hidden. Hide the graph.
    if (num_hidden === num_total) {
      plot_group_div.hide();
      return false;
    }
  }

  // Toggle buttons for y-axis limis
  // Handler for this is at top, so doesn't get created multiple times
  if (settings.ymax !== undefined || settings.ymin !== undefined) {
    let pgroup = $("#" + target).closest(".mqc_hcplot_plotgroup");
    let wrapper = $('<div class="mqc_hcplot_yaxis_limit_toggle hidden-xs" />').prependTo(pgroup);
    wrapper.append(
      '<span class="mqc_switch_wrapper" data-ymax="' +
        settings.ymax +
        '" data-ymin="' +
        settings.ymin +
        '" data-target="#' +
        target +
        '">Y-Limits: <span class="mqc_switch on">on</span></span>',
    );
    wrapper.after('<div class="clearfix" />');
  }

  // Translate this into JavaScript:
  // for sdata in view.data:
  //   if len(sdata["data"]) > 0 and isinstance(sdata["data"][0], list):
  //       x = [x[0] for x in sdata["data"]]
  //       y = [x[1] for x in sdata["data"]]
  //   else:
  //       x = [x for x in range(len(sdata["data"]))]
  //       y = sdata["data"]
  //
  //   fig.add_trace(
  //       go.Scatter(
  //           x=x,
  //           y=y,
  //           name=sdata["name"],
  //           mode="lines+markers",
  //           marker=dict(size=5),
  //       )
  //   )

  let traces = [];
  for (let sdata of data) {
    if (sdata.data.length > 0 && Array.isArray(sdata.data[0])) {
      x = sdata.data.map((x) => x[0]);
      y = sdata.data.map((x) => x[1]);
    } else {
      x = [...Array(sdata.data.length).keys()];
      y = sdata.data;
    }
    let trace = {
      type: "scatter",
      x: x,
      y: y,
      name: sdata.name,
      orientation: "h",
      marker: {
        line: {
          color: highlight_colors,
          width: highlight_colors.map((x) => (x ? 2 : 0)),
        },
      },
      marker_color: sdata.color,
    };
    traces.push(trace);
  }
  plot.datasets[dataset_idx].series = traces;
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
}
