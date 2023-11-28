////////////////////////////////////////////////
// Plotly Plotting Code
////////////////////////////////////////////////

// Global plot data variable. Accessed in many other JavaScript files.
let mqc_plots = {};

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

// Execute when page load has finished loading
$(function () {
  // Show loading warning
  let warning = $(".mqc_loading_warning").show();

  // Decompress the JSON plot data
  mqc_plots = JSON.parse(LZString.decompressFromBase64(mqc_compressed_plotdata));

  // Render plots on page load
  // (line below suppresses the "duplicated selector" linting inspection, as this time the
  // selector will return a different result)
  // noinspection JSJQueryEfficiency
  $(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)").each(function () {
    let target = $(this).attr("id");
    // Only one point per dataset, so multiply limit by arbitrary number.
    let max_num = mqc_config["num_datasets_plot_limit"] * 50;
    // Deferring each plot call prevents browser from locking up
    setTimeout(function () {
      plot_graph(target, undefined, max_num);
      if ($(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)").length === 0) {
        // All plots rendered successfully, so hiding the "loading" warning
        warning.hide();
      }
    }, 50);
  });

  // All plots already rendered successfully, so hiding the "loading" warning
  // noinspection JSJQueryEfficiency
  if ($(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)").length === 0) {
    warning.hide();
  }

  // Render a plot when clicked (heavy plots are not automatically rendered by default)
  let not_rendered = $(".hc-plot.not_rendered");
  $("body").on("click", ".render_plot", function (e) {
    let target = $(this).parent().attr("id");
    plot_graph(target);
    if (not_rendered.length === 0) {
      $("#mqc-warning-many-samples").hide();
    }
  });

  // Render all plots from header
  $("#mqc-render-all-plots").click(function () {
    not_rendered.each(function () {
      let target = $(this).attr("id");
      plot_graph(target);
    });
    $("#mqc-warning-many-samples").hide();
  });

  // Replot graphs when something changed in filters
  $(document).on("mqc_highlights mqc_renamesamples mqc_hidesamples", function () {
    // Replot graphs
    $(".hc-plot:not(.not_rendered)").each(function () {
      let target = $(this).attr("id");
      plot_graph(target);
    });
  });

  // A "Percentages" button above a plot is clicked
  $("button.switch_percent").click(function (e) {
    e.preventDefault();
    let target = $(this).data("target");

    // Toggling flags
    mqc_plots[target].p_active = !$(this).hasClass("active");
    $(this).toggleClass("active");

    // Replot graphs
    mqc_plots[target].replot();
    Plotly.relayout(target, "xaxis.tickformat", mqc_plots[target].p_active ? ".0%" : "");
  });

  // A "Log" button above a plot is clicked
  $("button.switch_log10").click(function (e) {
    e.preventDefault();
    let target = $(this).data("target");

    // Toggling flags
    mqc_plots[target].l_active = !$(this).hasClass("active");
    $(this).toggleClass("active");

    Plotly.relayout(target, "xaxis.type", mqc_plots[target].l_active ? "log" : "linear");
  });

  // Switch data source
  $(".dataset_switch_group button").click(function (e) {
    e.preventDefault();
    $(this).siblings("button.active").removeClass("active");
    $(this).addClass("active");
    let target = $(this).data("target");
    let active_dataset_idx = mqc_plots[target].active_dataset_idx;
    let new_dataset_idx = $(this).data("dataset_index");
    mqc_plots[target].active_dataset_idx = new_dataset_idx;
    if (active_dataset_idx === new_dataset_idx) return;
    mqc_plots[target].replot();
  });

  // Make divs height-draggable
  // http://jsfiddle.net/Lkwb86c8/
  $(".hc-plot:not(.no-handle)").each(function () {
    if (!$(this).parent().hasClass("hc-plot-wrapper")) {
      $(this).wrap('<div class="hc-plot-wrapper"></div>');
    }
    if (!$(this).siblings().hasClass("hc-plot-handle")) {
      $(this).after('<div class="hc-plot-handle"><span></span><span></span><span></span></div>');
    }
    $(this).css({ height: "auto", top: 0, bottom: "10px", position: "absolute" });
  });
  $(".hc-plot-handle").on("mousedown", function (e) {
    let wrapper = $(this).parent();
    let target = wrapper.children()[0].id;
    // $("#" + plot_id).find(".svg-container").height();
    // var plot_svg = $("#" + plot_id + " > .plot_container > .svg-container");
    // find .svg-container which is an immediate child of wrapper
    let startHeight = wrapper.height();
    let pY = e.pageY;
    $(document).on("mouseup", function () {
      // Clear listeners now that we've let go
      $(document).off("mousemove");
      $(document).off("mouseup");
      // Fire off a custom jQuery event for other javascript chunks to tie into
      // Bind to the plot div, which should have a custom ID
      $(wrapper.parent().find(".hc-plot, .beeswarm-plot")).trigger("mqc_plotresize");
    });
    $(document).on("mousemove", function (me) {
      var newHeight = startHeight + (me.pageY - pY);
      // container.css("height", newHeight);
      wrapper.css("height", newHeight);
      Plotly.relayout(target, { height: newHeight });
    });
  });

  // Special case for the beeswarm plot. TODO: fix for Plotly
  $(".hc-plot, .beeswarm-plot").on("mqc_plotresize", function (e) {
    if ($(this).highcharts()) {
      $(this).highcharts().reflow();
    }
  });

  // Switch the Y-axis limits on or off
  $(".mqc_hcplot_plotgroup").on("click", ".mqc_hcplot_yaxis_limit_toggle .mqc_switch_wrapper", function () {
    let target = $(this).data("target");
    let ymax = $(this).data("ymax");
    let ymin = $(this).data("ymin");
    let ymax_is_set = ymax !== "undefined" && ymax !== null;
    let ymin_is_set = ymin !== "undefined" && ymin !== null;
    if (!ymax_is_set && !(ymin_is_set && ymin !== 0))
      // If limits are not set (or the minimal limit is zero), don't do anything
      return;

    let y_limits_switch = $(this).find(".mqc_switch");
    if (y_limits_switch.hasClass("on")) {
      y_limits_switch.removeClass("on").addClass("off").text("off");
      Plotly.relayout(target, "yaxis.autorange", true);
      if (ymin === 0)
        // for plots with only positive numbers we want to keep the ymin=0 limit
        Plotly.relayout(target, "yaxis.range[0]", 0);
    } else {
      y_limits_switch.removeClass("off").addClass("on").text("on");
      Plotly.relayout(target, "yaxis.autorange", false);
      if (ymin_is_set && ymin !== 0) Plotly.relayout(target, "yaxis.range[0]", ymin);
      if (ymax_is_set) Plotly.relayout(target, "yaxis.range[1]", ymax);
    }
  });

  // Sort a heatmap by highlighted names
  $(".mqc_heatmap_sortHighlight").click(function (e) {
    e.preventDefault();
    let target = $(this).data("target").substr(1);
    if (mqc_plots[target]["config"]["sortHighlights"] === true) {
      mqc_plots[target]["config"]["sortHighlights"] = false;
      $(this).removeClass("active");
    } else {
      mqc_plots[target]["config"]["sortHighlights"] = true;
      $(this).addClass("active");
    }
    $(this).blur();
    plot_heatmap(target);
  });
});

// General function to rename samples. Takes a list of samples and a callback from
// plot-specific functions like plot_bar_graph() or plot_xy_line_graph()
function rename_samples(sample_names, rename_sample_func) {
  if (window.mqc_rename_f_texts.length > 0) {
    $.each(sample_names, function (sample_idx, sample_name) {
      $.each(window.mqc_rename_f_texts, function (idx, f_text) {
        let new_name;
        let t_text = window.mqc_rename_t_texts[idx];
        if (window.mqc_rename_regex_mode) {
          const re = new RegExp(f_text, "g");
          new_name = sample_name.replace(re, t_text);
        } else {
          new_name = sample_name.replace(f_text, t_text);
        }
        rename_sample_func(sample_idx, new_name);
      });
    });
  }
}

// General function to highlight samples
function get_highlight_colors(samples) {
  let highlight_colors = [];
  if (window.mqc_highlight_f_texts.length > 0) {
    $.each(samples, function (sample_idx, s_name) {
      highlight_colors[sample_idx] = null;
      $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
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
  return highlight_colors;
}

// General function to hide samples
function hide_samples(plot_group_div, samples, data) {
  plot_group_div.parent().find(".samples-hidden-warning").remove();
  plot_group_div.show();
  if (window.mqc_hide_f_texts.length > 0) {
    let num_hidden = 0;
    let num_total = samples.length;
    let j = samples.length;
    while (j--) {
      let match = false;
      for (let i = 0; i < window.mqc_hide_f_texts.length; i++) {
        const f_text = window.mqc_hide_f_texts[i];
        if (window.mqc_hide_regex_mode) {
          if (samples[j].match(f_text)) match = true;
        } else {
          if (samples[j].indexOf(f_text) > -1) match = true;
        }
      }
      if (window.mqc_hide_mode === "show") {
        match = !match;
      }
      if (match) {
        samples.splice(j, 1);
        $.each(data, function (k, d) {
          data.splice(j, 1);
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
}

// Call to render any plot
function plot_graph(target, dataset_idx, max_num) {
  let plot = mqc_plots[target];
  if (plot === undefined) return false;

  // If no dataset specified, check if there is an active button and set automatically
  if (dataset_idx === undefined) {
    let ds_buttons = $('.hc_switch_group button[data-action="set_data"][data-target="' + target + '"]');
    if (ds_buttons.length) dataset_idx = ds_buttons.filter(".active").data("dataset_index");
  }

  // If log status is not set and button is there, check whether it's active by default
  let config = plot["config"];
  if (config === undefined) {
    plot["config"] = {};
    config = {};
  }
  if (config["ytype"] === undefined) {
    const log_btn = $('.hc_switch_group button[data-action="set_log"][data-target="' + target + '"]');
    if (log_btn.length && log_btn.hasClass("active")) config["ytype"] = "logarithmic";
  }

  // XY Line charts
  if (plot.plot_type === "xy_line") {
    if (max_num === undefined || plot["datasets"][0].length < max_num) {
      plot_xy_line_graph(plot, target, dataset_idx);
      $("#" + target).removeClass("not_rendered");
    } else {
      $("#" + target)
        .addClass("not_rendered gt_max_num_ds")
        .html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
    }
  }
  // Bar graphs
  else if (plot.plot_type === "bar_graph") {
    if (max_num === undefined || plot["samples"][0].length < max_num) {
      plot_stacked_bar_graph(plot, target, dataset_idx);
      $("#" + target).removeClass("not_rendered");
    } else {
      $("#" + target)
        .addClass("not_rendered gt_max_num_ds")
        .html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
    }
  }
  // Scatter plots
  else if (plot["plot_type"] === "scatter") {
    if (max_num === undefined || Object.keys(plot["datasets"][0]).length < max_num) {
      plot_scatter_plot(plot, target, dataset_idx);
      $("#" + target).removeClass("not_rendered");
    } else {
      $("#" + target)
        .addClass("not_rendered gt_max_num_ds")
        .html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
    }
  }
  // Beeswarm graphs
  else if (plot["plot_type"] === "beeswarm") {
    if (max_num === undefined || plot["samples"][0].length < max_num) {
      plot_beeswarm_graph(plot, target, dataset_idx);
      $("#" + target).removeClass("not_rendered");
    } else {
      $("#" + target)
        .addClass("not_rendered gt_max_num_ds")
        .html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
    }
  }
  // Heatmap plots
  else if (plot["plot_type"] === "heatmap") {
    if (max_num === undefined || plot["xcats"][0].length < max_num) {
      plot_heatmap(plot, target, dataset_idx);
      $("#" + target).removeClass("not_rendered");
    } else {
      $("#" + target)
        .addClass("not_rendered gt_max_num_ds")
        .html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
    }
  }
  // Not recognised
  else {
    console.log("Did not recognise plot type: " + plot["plot_type"]);
  }
}

// Highlight text with a fadeout background colour highlight
function highlight_fade_text(obj) {
  let orig_col = $(obj).css("color");
  obj.css({
    display: "inline-block",
    "background-color": "#5bc0de",
    color: "#FFFFFF",
    WebkitTransition: "background-color 0s, color 0s",
    MozTransition: "background-color 0s, color 0s",
    MsTransition: "background-color 0s, color 0s",
    OTransition: "background-color 0s, color 0s",
    transition: "background-color 0s, color 0s",
  });
  setTimeout(function () {
    obj.css({
      "background-color": "#FFFFFF",
      color: orig_col,
      WebkitTransition: "background-color 0.5s, color 0.5s",
      MozTransition: "background-color 0.5s, color 0.5s",
      MsTransition: "background-color 0.5s, color 0.5s",
      OTransition: "background-color 0.5s, color 0.5s",
      transition: "background-color 0.5s, color 0.5s",
    });
  }, 500);
}
