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
    this.target = dump["id"];
    this.datasets = dump["datasets"];
    this.square = dump["square"];
    this.axisControlledBySwitches = dump["axis_controlled_by_switches"];
    this.static = dump["static"] ?? false;
    // To make sure we only render plot once
    this.rendered = false;
    // State of toggles
    this.activeDatasetIdx = 0;
    this.lActive = dump["l_active"];
    this.pActive = dump["p_active"];
  }

  activeDatasetSize() {
    throw new Error("activeDatasetSize() not implemented");
  }

  resize(newHeight) {
    if (!this.square) {
      Plotly.relayout(this.target, { autosize: true, height: newHeight });
    } else {
      // noinspection JSSuspiciousNameCombination
      Plotly.relayout(this.target, { autosize: true, height: newHeight, width: newHeight });
    }
  }

  buildTraces(layout) {
    throw new Error("buildTraces() not implemented");
  }

  afterPlotCreated() {
    // Do nothing
  }
}

function initPlot(dump) {
  if (dump["plot_type"] === "xy_line") return new LinePlot(dump);
  if (dump["plot_type"] === "bar_graph") return new BarPlot(dump);
  if (dump["plot_type"] === "scatter") return new ScatterPlot(dump);
  if (dump["plot_type"] === "heatmap") return new HeatmapPlot(dump);
  if (dump["plot_type"] === "violin") return new ViolinPlot(dump);
  console.log("Did not recognise plot type: " + dump["plot_type"]);
  return null;
}

// Execute when page load has finished loading
$(function () {
  // Show loading warning
  let loading_warning = $(".mqc_loading_warning").show();

  // Decompress the JSON plot data and init plot objects
  let mqc_plotdata = JSON.parse(LZString.decompressFromBase64(mqc_compressed_plotdata));
  mqc_plots = Object.fromEntries(Object.values(mqc_plotdata).map((data) => [data.id, initPlot(data)]));

  let shouldRender = $(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)");

  // Render plots on page load
  shouldRender.each(function () {
    let target = $(this).attr("id");
    // Deferring each plot call prevents browser from locking up
    setTimeout(function () {
      renderPlot(target);
      if ($(".hc-plot.not_rendered:visible:not(.gt_max_num_ds)").length === 0)
        // All plots rendered successfully (or hidden with gt_max_num_ds), so hiding the warning
        $(".mqc_loading_warning").hide();
    }, 50);
  });

  // All plots rendered successfully (or hidden with gt_max_num_ds), so hiding the warning
  if (shouldRender.length === 0) loading_warning.hide();

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
    mqc_plots[target].pActive = !$(this).hasClass("active");
    $(this).toggleClass("active");

    renderPlot(target);
  });

  // A "Log" button above a plot is clicked
  $("button.interactive-switch-group.log10-switch").click(function (e) {
    e.preventDefault();
    let target = $(this).data("pid");

    // Toggling flags
    mqc_plots[target].lActive = !$(this).hasClass("active");
    $(this).toggleClass("active");

    mqc_plots[target].axisControlledBySwitches.map((axis) => {
      Plotly.relayout(target, axis + ".type", mqc_plots[target].lActive ? "log" : "linear");
    });
  });

  // Switch data source
  $(".interactive-switch-group.dataset-switch-group button").click(function (e) {
    e.preventDefault();
    if ($(this).hasClass("active")) return;
    $(this).siblings("button.active").removeClass("active");
    $(this).addClass("active");
    let target = $(this).data("pid");
    let activeDatasetIdx = mqc_plots[target].activeDatasetIdx;
    let newDatasetIdx = $(this).data("datasetIndex");
    mqc_plots[target].activeDatasetIdx = newDatasetIdx;
    if (activeDatasetIdx === newDatasetIdx) return;

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
    $(this).css({ height: "auto", top: 0, bottom: "6px", position: "absolute" });
  });

  $(".hc-plot-handle").on("mousedown", function (e) {
    let wrapper = $(this).parent();
    let target = wrapper.children(".hc-plot")[0].id;
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
      let newHeight = startHeight + (me.pageY - pY) + 2; // 2 px for the border or something
      wrapper.css("height", newHeight);
      if (mqc_plots[target] !== undefined) mqc_plots[target].resize(newHeight - 7); // 7 is the height of the handle overlapping the plot wrapper
    });
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

  // Useful for default views like a table that can be switched to a Violin plot.
  if (plot.static) return;

  let container = $("#" + target);

  // When the plot was already rendered, it's faster to call react, using the same signature
  // https://plotly.com/javascript/plotlyjs-function-reference/#plotlyreact
  let func;
  if (!plot.rendered) {
    func = Plotly.newPlot;
    plot.rendered = true;
    container.removeClass("not_rendered").parent().find(".render_plot").remove();
    if ($(".hc-plot.not_rendered").length === 0) $("#mqc-warning-many-samples").hide();
  } else {
    func = Plotly.react;
  }

  let dataset = plot.datasets[plot.activeDatasetIdx];
  let layout = JSON.parse(JSON.stringify(dataset.layout));

  // Apply pct/log toggle states
  plot.axisControlledBySwitches.map((axis) => {
    layout[axis].type = plot.lActive ? "log" : "linear";
  });
  if (plot.pActive) {
    layout.xaxis.ticksuffix = "%";
    layout.xaxis.hoverformat = ".1f";
    layout.xaxis.range = [0, 100];
  }

  let traces = plot.buildTraces(layout);
  if (traces.length > 0 && traces[0].constructor === Array) traces = [].concat.apply([], traces); // if list of lists, flatten
  if (traces.length === 0) {
    // All series hidden. Hide the graph.
    container.hide();
    return;
  }

  container.show();
  func(target, traces, layout, {
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
      "toImage",
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

function updateObject(target, source) {
  // Iterate through all keys in the source object
  for (const key in source) {
    // Check if the value is not null
    if (source[key] !== null) {
      // Check if the value is an object and not an array
      if (typeof source[key] === "object" && !Array.isArray(source[key])) {
        // If the target doesn't have this key or it's not an object, initialize it
        if (!target[key] || typeof target[key] !== "object") {
          target[key] = {};
        }
        // Recursively update the object
        updateObject(target[key], source[key]);
      } else {
        // Directly update the value
        target[key] = source[key];
      }
    }
  }
}