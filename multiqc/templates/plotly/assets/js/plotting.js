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
window.mqc_hide_mode = "hide";
window.mqc_hide_f_texts = [];
window.mqc_hide_regex_mode = false;

class Plot {
  constructor(dump) {
    this.target = dump.id;
    this.plot_type = dump.plot_type;
    this.layout = dump.layout;
    this.trace_params = dump.trace_params;
    this.datasets = dump.datasets;
    this.config = dump.config;
    this.axis_controlled_by_switches = dump.axis_controlled_by_switches;
    // To make sure we only render plot once
    this.rendered = false;
    // State of toggles
    this.active_dataset_idx = 0;
    this.p_active = dump.config.p_active;
    this.l_active = dump.config.l_active;
  }

  activeDatasetSize() {
    throw new Error("activeDatasetSize() not implemented");
  }

  resize(newHeight) {
    Plotly.relayout(this.target, { height: newHeight });
    if (this.config.square) {
      Plotly.relayout(this.target, { width: newHeight });
    }
  }

  buildTraces() {
    throw new Error("buildTraces() not implemented");
  }

  afterPlotCreated() {
    // Do nothing
  }
}

function initPlot(dump) {
  if (dump.plot_type === "xy_line") return new LinePlot(dump);
  if (dump.plot_type === "bar_graph") return new BarPlot(dump);
  if (dump.plot_type === "scatter") return new ScatterPlot(dump);
  if (dump.plot_type === "heatmap") return new HeatmapPlot(dump);
  if (dump.plot_type === "violin") return new ViolinPlot(dump);
  console.log("Did not recognise plot type: " + dump.plot_type);
  return null;
}

// Execute when page load has finished loading
$(function () {
  // Show loading warning
  let loading_warning = $(".mqc_loading_warning").show();

  // Decompress the JSON plot data and init plot objects
  let mqc_plotdata = JSON.parse(LZString.decompressFromBase64(mqc_compressed_plotdata));
  mqc_plots = Object.fromEntries(Object.values(mqc_plotdata).map((data) => [data.id, initPlot(data)]));

  let should_render = $(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)");

  // Render plots on page load
  should_render.each(function () {
    let target = $(this).attr("id");
    // Only one point per dataset, so multiply limit by arbitrary number.
    let maxNum = mqc_config["num_datasets_plot_limit"] * 50;
    // Deferring each plot call prevents browser from locking up
    setTimeout(function () {
      let plot = mqc_plots[target];
      if (plot.activeDatasetSize() > maxNum) {
        $("#" + target)
          .addClass("not_rendered gt_max_num_ds")
          .html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
      } else {
        renderPlot(target);
        if ($(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)").length === 0)
          // All plots rendered successfully (or hidden with gt_max_num_ds), so hiding the warning
          $(".mqc_loading_warning").hide();
      }
    }, 50);
  });

  // All plots rendered successfully (or hidden with gt_max_num_ds), so hiding the warning
  if (should_render.length === 0) loading_warning.hide();

  // Render a plot when clicked (heavy plots are not automatically rendered by default)
  $("body").on("click", ".render_plot", function (e) {
    renderPlot($(this).parent().attr("id"));
  });

  // Render all plots from header, even those that are hidden
  $("#mqc-render-all-plots").click(function () {
    $(".hc-plot.not_rendered").each(function () {
      renderPlot($(this).attr("id"));
    });
  });

  // Replot graphs when something changed in filters
  $(document).on("mqc_highlights mqc_renamesamples mqc_hidesamples", function () {
    // Replot graphs
    $(".hc-plot:not(.not_rendered)").each(function () {
      renderPlot($(this).attr("id"));
    });
  });

  // A "Percentages" button above a plot is clicked
  $("button.interactive-switch-group.percent-switch").click(function (e) {
    e.preventDefault();
    let target = $(this).data("pid");

    // Toggling flags
    mqc_plots[target].p_active = !$(this).hasClass("active");
    $(this).toggleClass("active");

    renderPlot(target);
  });

  // A "Log" button above a plot is clicked
  $("button.interactive-switch-group.log10-switch").click(function (e) {
    e.preventDefault();
    let target = $(this).data("pid");

    // Toggling flags
    mqc_plots[target].l_active = !$(this).hasClass("active");
    $(this).toggleClass("active");

    mqc_plots[target].axis_controlled_by_switches.map((axis) => {
      Plotly.relayout(target, axis + ".type", mqc_plots[target].l_active ? "log" : "linear");
    });
  });

  // Switch data source
  $(".interactive-switch-group.dataset-switch-group button").click(function (e) {
    e.preventDefault();
    if ($(this).hasClass("active")) return;
    $(this).siblings("button.active").removeClass("active");
    $(this).addClass("active");
    let target = $(this).data("pid");
    let active_dataset_idx = mqc_plots[target].active_dataset_idx;
    let new_dataset_idx = $(this).data("datasetIndex");
    mqc_plots[target].active_dataset_idx = new_dataset_idx;
    if (active_dataset_idx === new_dataset_idx) return;

    renderPlot(target);
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
      let newHeight = startHeight + (me.pageY - pY);
      wrapper.css("height", newHeight);
      mqc_plots[target].resize(newHeight);
    });
  });

  // // Special case for the beeswarm plot. TODO: fix for Plotly
  // $(".hc-plot, .beeswarm-plot").on("mqc_plotresize", function (e) {
  //   if ($(this).highcharts()) {
  //     $(this).highcharts().reflow();
  //   }
  // });

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

  // Sort a heatmap by highlighted names  // TODO: fix for Plotly
  $(".mqc_heatmap_sortHighlight").click(function (e) {
    e.preventDefault();
    let target = $(this).data("target").substr(1);
    if (mqc_plots[target].sort_highlights === true) {
      mqc_plots[target].sort_highlights = false;
      $(this).removeClass("active");
    } else {
      mqc_plots[target].sort_highlights = true;
      $(this).addClass("active");
    }
    $(this).blur();
    renderPlot(target);
  });
});

// Highlighting, hiding and renaming samples. Takes a list of samples, returns
// a list of objects: {"name": "new_name", "highlight": "#cccccc", "hidden": false}
function applyToolboxSettings(samples, target) {
  // init object with default values
  let objects = samples.map((name) => ({ name: name, highlight: null, hidden: false }));

  // Rename samples
  if (window.mqc_rename_f_texts.length > 0) {
    objects.map((obj) => {
      for (let p_idx = 0; p_idx < window.mqc_rename_f_texts.length; p_idx++) {
        let pattern = window.mqc_rename_f_texts[p_idx];
        let new_text = window.mqc_rename_t_texts[p_idx];
        obj.name = obj.name.replace(pattern, new_text);
      }
    });
  }

  // Highlight samples
  if (window.mqc_highlight_f_texts.length > 0) {
    objects.map((obj) => {
      for (let i = 0; i < window.mqc_highlight_f_texts.length; i++) {
        const f_text = window.mqc_highlight_f_texts[i];
        const f_col = window.mqc_highlight_f_cols[i];
        let match = false;
        if (window.mqc_highlight_regex_mode) {
          if (obj.name.match(f_text)) match = true;
        } else {
          if (obj.name.indexOf(f_text) > -1) match = true;
        }
        if (match) obj.highlight = f_col;
      }
    });
  }

  // Hide samples
  if (window.mqc_hide_f_texts.length > 0) {
    let groupDiv = $("#" + target).closest(".mqc_hcplot_plotgroup");
    groupDiv.parent().find(".samples-hidden-warning").remove();
    groupDiv.show();

    objects.map((obj) => {
      let match = false;
      for (let i = 0; i < window.mqc_hide_f_texts.length; i++) {
        const f_text = window.mqc_hide_f_texts[i];
        if (window.mqc_hide_regex_mode) {
          if (obj.name.match(f_text)) match = true;
        } else {
          if (obj.name.indexOf(f_text) > -1) match = true;
        }
      }
      if (window.mqc_hide_mode === "show") match = !match;
      if (match) obj.hidden = true;
    });

    // Some series hidden. Show a warning text string.
    let nHidden = objects.filter((obj) => obj.hidden).length;
    if (nHidden > 0) {
      const alert =
        '<div class="samples-hidden-warning alert alert-warning">' +
        '<span class="glyphicon glyphicon-info-sign"></span>' +
        "<strong>Warning:</strong> " +
        nHidden +
        " samples hidden. " +
        '<a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a>' +
        "</div>";
      groupDiv.before(alert);
    }
    // All series hidden. Hide the graph.
    if (nHidden === objects.length) {
      groupDiv.hide();
      return null;
    }
  }

  // Return the object indexed by sample names
  return objects;
}

// // Hiding samples. Returns indices of samples in the "samples" array
// function hideSamples(target, samples) {
//   let plot_group_div = $("#" + target).closest(".mqc_hcplot_plotgroup");
//   plot_group_div.parent().find(".samples-hidden-warning").remove();
//   plot_group_div.show();
//
//   if (window.mqc_hide_f_texts.length === 0) return [];
//
//   let result = [];
//   for (let j = 0; j < samples.length; j++) {
//     let match = false;
//     for (let i = 0; i < window.mqc_hide_f_texts.length; i++) {
//       const f_text = window.mqc_hide_f_texts[i];
//       if (window.mqc_hide_regex_mode) {
//         if (samples[j].match(f_text)) match = true;
//       } else {
//         if (samples[j].indexOf(f_text) > -1) match = true;
//       }
//     }
//     if (window.mqc_hide_mode === "show") {
//       match = !match;
//     }
//     if (match) {
//       result.push(j);
//     }
//   }
//   // Some series hidden. Show a warning text string.
//   if (result.length > 0) {
//     const alert =
//       '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> ' +
//       result.length +
//       ' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
//     plot_group_div.before(alert);
//   }
//   // All series hidden. Hide the graph.
//   if (result.length === samples.length) plot_group_div.hide();
//   return result;
// }

// Call to render any plot
function renderPlot(target) {
  let plot = mqc_plots[target];
  if (plot === undefined) return false;
  if (plot.datasets.length === 0) return false;

  // When the plot was already rendered, it's faster to call react, using the same signature
  // https://plotly.com/javascript/plotlyjs-function-reference/#plotlyreact
  let func;
  if (!plot.rendered) {
    func = Plotly.newPlot;
    plot.rendered = true;
    $("#" + target)
      .removeClass("not_rendered")
      .parent()
      .find(".render_plot")
      .remove();
    if ($(".hc-plot.not_rendered").length === 0) $("#mqc-warning-many-samples").hide();
  } else {
    func = Plotly.react;
  }

  // Apply toggle states
  plot.layout.xaxis.tickformat = plot.p_active ? ".0%" : "";
  plot.axis_controlled_by_switches.map((axis) => {
    plot.layout[axis].type = plot.l_active ? "log" : "linear";
  });

  func(target, plot.buildTraces(), plot.layout, {
    responsive: true,
    displaylogo: false,
    displayModeBar: true,
    toImageButtonOptions: { filename: target },
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

  plot.afterPlotCreated();
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
