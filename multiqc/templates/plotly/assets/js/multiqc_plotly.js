////////////////////////////////////////////////
// Plotly Plotting Code
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

// Execute when page load has finished loading
$(function () {
  // Show loading warning
  $(".mqc_loading_warning").show();

  // Decompress the JSON plot data
  mqc_plots = JSON.parse(LZString.decompressFromBase64(mqc_compressed_plotdata));

  // Render plots on page load
  $(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)").each(function () {
    var target = $(this).attr("id");
    // Only one point per dataset, so multiply limit by arbitrary number.
    var max_num = mqc_config["num_datasets_plot_limit"] * 50;
    // Deferring each plot call prevents browser from locking up
    setTimeout(function () {
      plot_graph(target, undefined, max_num);
      if ($(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)").length == 0) {
        $(".mqc_loading_warning").hide();
      }
    }, 50);
  });
  if ($(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)").length == 0) {
    $(".mqc_loading_warning").hide();
  }

  // Render a plot when clicked
  $("body").on("click", ".render_plot", function (e) {
    var target = $(this).parent().attr("id");
    plot_graph(target);
    if ($(".hc-plot.not_rendered").length == 0) {
      $("#mqc-warning-many-samples").hide();
    }
  });

  // Render all plots from header
  $("#mqc-render-all-plots").click(function () {
    $(".hc-plot.not_rendered").each(function () {
      var target = $(this).attr("id");
      plot_graph(target);
    });
    $("#mqc-warning-many-samples").hide();
  });

  // Replot graphs when something changed in filters
  $(document).on("mqc_highlights mqc_renamesamples mqc_hidesamples", function () {
    // Replot graphs
    $(".hc-plot:not(.not_rendered)").each(function () {
      var target = $(this).attr("id");
      plot_graph(target);
    });
  });

  $("button.switch_percent").click(function (e) {
    e.preventDefault();
    let target = $(this).data("target");
    let active_dataset_idx = mqc_plots[target].active_dataset_idx;
    let dataset = mqc_plots[target].datasets[active_dataset_idx];

    // Toggling flags
    mqc_plots[target].p_active = !$(this).hasClass("active");
    $(this).toggleClass("active");

    let x = [];
    for (let cat of dataset) {
      x.push(mqc_plots[target].p_active ? cat.data_pct : cat.data);
    }
    Plotly.restyle(target, "x", x);
    Plotly.relayout(target, "xaxis.tickformat", mqc_plots[target].p_active ? ".0%" : "");
  });

  $("button.switch_log10").click(function (e) {
    e.preventDefault();
    var target = $(this).data("target");

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
    var target = $(this).data("target");
    var active_dataset_idx = mqc_plots[target].active_dataset_idx;
    var new_dataset_idx = $(this).data("dataset_index");
    mqc_plots[target].active_dataset_idx = new_dataset_idx;
    if (active_dataset_idx === new_dataset_idx) {
      return;
    }
    var dataset = mqc_plots[target]["datasets"][new_dataset_idx];

    let x = [];
    for (let cat of dataset) {
      x.push(dataset.p_active ? cat.data_pct : cat.data);
    }
    Plotly.restyle(target, "x", x);

    const ymax = $(this).data("xmax");
    if (ymax) {
      Plotly.relayout(target, "yaxis.range", [null, ymax]);
    }
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
    var wrapper = $(this).parent();
    var target = wrapper.children()[0].id;
    // $("#" + plot_id).find(".svg-container").height();
    // var plot_svg = $("#" + plot_id + " > .plot_container > .svg-container");
    // find .svg-container which is an immediate child of wrapper
    var startHeight = wrapper.height();
    var pY = e.pageY;
    $(document).on("mouseup", function (e) {
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

  $(".hc-plot, .beeswarm-plot").on("mqc_plotresize", function (e) {
    if ($(this).highcharts()) {
      $(this).highcharts().reflow();
    }
  });

  // Switch a y-axis limit on or off
  $(".mqc_hcplot_plotgroup").on("click", ".mqc_hcplot_yaxis_limit_toggle .mqc_switch_wrapper", function () {
    var target = $($(this).data("target")).highcharts();
    var ymax = $(this).data("ymax");
    var ymin = $(this).data("ymin");
    ymax = ymax == "undefined" ? null : ymax;
    ymin = ymin == "undefined" ? null : ymin;
    var mqc_switch = $(this).find(".mqc_switch");
    if (mqc_switch.hasClass("on")) {
      target.yAxis[0].update({ max: null, min: null });
      mqc_switch.removeClass("on").addClass("off").text("off");
    } else {
      target.yAxis[0].update({ max: ymax, min: ymin });
      mqc_switch.removeClass("off").addClass("on").text("on");
    }
  });

  // Sort a heatmap by highlighted names
  $(".mqc_heatmap_sortHighlight").click(function (e) {
    e.preventDefault();
    var target = $(this).data("target").substr(1);
    if (mqc_plots[target]["config"]["sortHighlights"] == true) {
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

// Call to render any plot
function plot_graph(target, dataset_idx, max_num) {
  if (mqc_plots[target] === undefined) {
    return false;
  }
  let plot = mqc_plots[target];

  // If no dataset specified, check if there is an active button and set automatically
  if (dataset_idx === undefined) {
    var ds_btns = $('.hc_switch_group button[data-action="set_data"][data-target="' + target + '"]');
    if (ds_btns.length) {
      dataset_idx = ds_btns.filter(".active").data("dataset_index");
    }
  }

  // If log status is not set and button is there, check whether it's active by default
  let config = plot["config"];
  if (config === undefined) {
    plot["config"] = {};
    config = {};
  }
  if (config["ytype"] === undefined) {
    const log_btn = $('.hc_switch_group button[data-action="set_log"][data-target="' + target + '"]');
    if (log_btn.length && log_btn.hasClass("active")) {
      config["ytype"] = "logarithmic";
    }
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
  var orig_col = $(obj).css("color");
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
