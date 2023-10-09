// Stacked Bar Graph
function plot_stacked_bar_graph(plot, target, dataset_idx) {
  if (plot === undefined || plot["plot_type"] !== "bar_graph") {
    return false;
  }
  if (dataset_idx === undefined) {
    dataset_idx = 0;
  }

  // Make a clone of everything, so that we can mess with it,
  // while keeping the original data intact
  let data = JSON.parse(JSON.stringify(plot["datasets"][dataset_idx]));
  let samples = JSON.parse(JSON.stringify(plot["samples"][dataset_idx]));
  let layout = JSON.parse(JSON.stringify(plot["layout"]));

  // Rename samples
  if (window.mqc_rename_f_texts.length > 0) {
    $.each(samples, function (j, s_name) {
      $.each(window.mqc_rename_f_texts, function (idx, f_text) {
        if (window.mqc_rename_regex_mode) {
          const re = new RegExp(f_text, "g");
          samples[j] = samples[j].replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          samples[j] = samples[j].replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
    });
  }

  // Highlight samples
  let highlight_colors = [];
  if (window.mqc_highlight_f_texts.length > 0) {
    $.each(samples, function (sample_idx, s_name) {
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
    let num_total = samples.length;
    let j = samples.length;
    while (j--) {
      let s_name = samples[j];
      let match = false;
      for (let i = 0; i < window.mqc_hide_f_texts.length; i++) {
        const f_text = window.mqc_hide_f_texts[i];
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
      if (window.mqc_hide_mode === "show") {
        match = !match;
      }
      if (match) {
        samples.splice(j, 1);
        $.each(data, function (k, d) {
          data[k]["data"].splice(j, 1);
        });
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

  // Render the plotly plot
  let series = [];
  for (let cat of data) {
    let trace = {
      type: "bar",
      y: samples,
      x: cat.data,
      name: cat.name,
      orientation: "h",
      marker: {
        line: {
          color: highlight_colors,
          width: highlight_colors.map((x) => (x ? 2 : 0)),
        },
      },
      marker_color: cat.color,
    };
    series.push(trace);
  }
  plot.datasets[dataset_idx].series = series;
  Plotly.newPlot(target, series, layout, {
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
