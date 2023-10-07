////////////////////////////////////////////////
// HighCharts Plotting Code
////////////////////////////////////////////////

// Global plot data variable
mqc_plots = {};

// Initialise the toolbox filters
window.mqc_highlight_f_texts = [];
window.mqc_highlight_f_cols = [];
window.mqc_highlight_regex_mode = false;
window.mqc_rename_f_texts = [];
window.mqc_rename_t_texts = [];
window.mqc_rename_regex_mode = false;
window.mqc_hide_mode = "hide";
window.mqc_hide_f_texts = [];
window.mqc_hide_regex_mode = false;
window.HCDefaults = undefined;

// Stacked Bar Graph
function plot_stacked_bar_graph(target, ds) {
  if (mqc_plots[target] === undefined || mqc_plots[target]["plot_type"] !== "bar_graph") {
    return false;
  }
  if (ds === undefined) {
    ds = 0;
  }

  // Make a clone of everything, so that we can mess with it,
  // while keeping the original data intact
  var data = JSON.parse(JSON.stringify(mqc_plots[target]["datasets"][ds]));
  var samples = JSON.parse(JSON.stringify(mqc_plots[target]["samples"][ds]));
  var config = JSON.parse(JSON.stringify(mqc_plots[target]["config"]));
  var layout = JSON.parse(JSON.stringify(mqc_plots[target]["layout"]));

  if (config["stacking"] === undefined) {
    config["stacking"] = "normal";
  }
  if (config["stacking"] === "normal") {
    config["groupPadding"] = "0.02";
    config["pointPadding"] = "0.1";
  } else {
    config["groupPadding"] = "0.1";
    config["pointPadding"] = 0;
  }
  if (config["ytype"] === undefined) {
    config["ytype"] = "linear";
  }
  if (config["reversedStacks"] === undefined) {
    config["reversedStacks"] = false;
  }
  if (config["use_legend"] === undefined) {
    config["use_legend"] = true;
  }
  if (config["yDecimals"] === undefined) {
    config["yDecimals"] = true;
  }
  if (config["click_func"] === undefined) {
    config["click_func"] = function () {};
  } else {
    if (config["cursor"] === undefined) {
      config["cursor"] = "pointer";
    }
  }
  if (config["tt_percentages"] === undefined) {
    config["tt_percentages"] = true;
  }
  if (config["borderWidth"] === undefined) {
    config["borderWidth"] = 0;
  }

  if (config["ytype"] == "logarithmic") {
    if (config["ymin"] == 0 || config["ymin"] == undefined) {
      config["ymin"] = 1;
    }
    var minTickInt = "auto";
  } else {
    var minTickInt = undefined;
  }

  // Rename samples
  if (window.mqc_rename_f_texts.length > 0) {
    $.each(samples, function (j, s_name) {
      $.each(window.mqc_rename_f_texts, function (idx, f_text) {
        if (window.mqc_rename_regex_mode) {
          var re = new RegExp(f_text, "g");
          samples[j] = samples[j].replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          samples[j] = samples[j].replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
    });
  }

  // Highlight samples
  if (window.mqc_highlight_f_texts.length > 0) {
    $.each(samples, function (j, s_name) {
      $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
        if (f_text == "") {
          return true;
        } // skip blanks
        if (
          (window.mqc_highlight_regex_mode && s_name.match(f_text)) ||
          (!window.mqc_highlight_regex_mode && s_name.indexOf(f_text) > -1)
        ) {
          // Make the data point in each series with this index have a border colour
          $.each(data, function (k, d) {
            data[k]["data"][j] = {
              y: data[k]["data"][j],
              borderColor: window.mqc_highlight_f_cols[idx],
            };
          });
        }
      });
    });
    // Bump the borderWidth to make the highlights more obvious
    if (config["borderWidth"] <= 2) {
      config["borderWidth"] = 2;
    }
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
    var num_total = samples.length;
    var j = samples.length;
    while (j--) {
      var s_name = samples[j];
      var match = false;
      for (i = 0; i < window.mqc_hide_f_texts.length; i++) {
        var f_text = window.mqc_hide_f_texts[i];
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
        samples.splice(j, 1);
        $.each(data, function (k, d) {
          data[k]["data"].splice(j, 1);
        });
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

  // Make the plotly plot
  let series = [];
  for (let cat of data) {
    let trace = {
      type: "bar",
      y: samples,
      x: cat.data,
      name: cat.name,
      orientation: "h",
      marker: {
        color: cat.color,
        line: { width: 0 },
      },
    };
    series.push(trace);
  }

  mqc_plots[target].datasets[ds].series = series;
  Plotly.newPlot(target, series, layout);
}
