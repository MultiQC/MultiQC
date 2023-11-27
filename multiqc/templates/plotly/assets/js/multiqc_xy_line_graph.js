// Basic Line Graph
function plot_xy_line_graph(plot, target, dataset_idx) {
  if (plot === undefined || plot["plot_type"] !== "xy_line") return false;

  if (dataset_idx === undefined) dataset_idx = 0;

  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data intact
  let data = JSON.parse(JSON.stringify(plot["datasets"][dataset_idx]));
  let layout = JSON.parse(JSON.stringify(plot["layout"]));
  let pconfig = JSON.parse(JSON.stringify(plot["pconfig"]));

  // Rename samples
  let _s_names = [];
  for (let sdata of data) {
    _s_names.push(sdata.name);
  }
  rename_samples(_s_names, function (sample_idx, new_name) {
    data[sample_idx]["name"] = new_name;
  });

  // Highlight samples
  let highlight_colors = [];
  if (window.mqc_highlight_f_texts.length > 0) {
    $.each(data, function (sample_idx) {
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
  for (let sdata of data) {
    let [x, y] = sdata_to_xy(sdata);
    let trace = {
      type: "scatter",
      x: x,
      y: y,
      name: sdata.name,
      orientation: "h",
      mode: "lines",
      marker: {
        size: 5,
        line: {
          color: highlight_colors,
          width: highlight_colors.map((x) => (x ? 2 : 0)),
        },
      },
      marker_color: sdata.color,
      // hovertemplate: pconfig["tt_label"]
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

  // Function to re-render plot with new data
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

    // No need to restyle because the data is already limited to requested xmax and ymax in Python
    // if ($(this).data("xmax") || $(this).data("ymax")) {
    //   Plotly.relayout(target, "yaxis.range", [$(this).data("xmax") || null, $(this).data("ymax") || null]);
    // }
  };
}
