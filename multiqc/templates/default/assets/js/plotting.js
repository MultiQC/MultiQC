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
    this.anchor = dump["anchor"];
    this.layout = dump["layout"];
    this.datasets = dump["datasets"];
    this.pctAxisUpdate = dump["pct_axis_update"];
    this.axisControlledBySwitches = dump["axis_controlled_by_switches"];
    this.square = dump["square"];
    // To make sure we only render plot once
    this.rendered = false;
    // State of toggles
    this.activeDatasetIdx = 0;
    this.lActive = dump["l_active"];
    this.pActive = dump["p_active"];
    this.deferRender = dump["defer_render"];
  }

  activeDatasetSize() {
    throw new Error("activeDatasetSize() not implemented");
  }

  resize(newHeight, newWidth) {
    if (newHeight === null || newHeight === undefined) {
      console.error("Plot.resize: newHeight is " + newHeight);
      return;
    }
    this.layout.height = newHeight;
    this.layout.autosize = true;

    if (newWidth !== null && newWidth !== undefined) {
      this.layout.width = newWidth;
      this.layout.autosize = false;
    } else if (this.square) {
      // noinspection JSSuspiciousNameCombination
      this.layout.width = newHeight;
    } else {
      this.layout.width = null;
    }
    Plotly.relayout(this.anchor, this.layout);
  }

  buildTraces() {
    throw new Error("buildTraces() not implemented");
  }

  afterPlotCreated() {
    // Do nothing
  }

  recalculateTicks(filteredSettings, axis, maxTicks) {
    // Recalculate tick settings (array or auto; highlight color) for the axis,
    // based on the desired number of ticks and the sample toolbox settings.
    function subsample(values, num, start = 0, roundBin = true) {
      // Take ~`num` samples from values evenly, always include `start`.
      // If `roundBin` is true, the bins will be rounded to the nearest integer,
      // so the ticks will be evenly distributed, but the total number of ticks
      // may be less than `num`.

      if (values.length <= num) return values;
      if (values.length <= 1) return values;
      if (num === 0) return [];
      if (num === 1) return [values[start]];

      let binSize = (values.length - 1) / (num - 1);
      if (roundBin) binSize = Math.ceil(binSize);

      // Split into two halves: before and after pivot, including pivot into both. This way
      // we want to make sure pivot is always included in the result.
      let indices = Array.from({ length: values.length }, (_, i) => i);
      let after = indices.slice(start);
      let before = indices.slice(0, start + 1); // including the pivot

      // Stepping forward `binsize` steps, starting from the pivot
      after = Array.from({ length: after.length }, (_, i) => Math.ceil(binSize * i))
        .filter((index) => index < after.length)
        .map((index) => after[index]);

      before.reverse(); // Stepping back starting from the pivot
      before = Array.from({ length: before.length }, (_, i) => Math.ceil(binSize * i))
        .filter((index) => index < before.length)
        .map((index) => before[index]);
      before.reverse();
      before = before.slice(0, before.length - 1); // remove the pivot

      indices = before.concat(after);
      return indices.map((i) => values[i]);
    }

    let highlighted = filteredSettings.filter((s) => s.highlight);
    let firstHighlightedSample = this.firstHighlightedSample(filteredSettings);

    if (highlighted.length === 0) {
      axis.tickmode = null;
      axis.tickvals = null;
      axis.ticktext = null;
    } else {
      // Have to switch to tickmode=array to set colors to ticks. however, this way plotly will try
      // to fit _all_ ticks on the screen, and if there are too many, they will overlap. to prevent that,
      // if there are too many samples, we will show only highlighted samples plus a subsampled number
      // of ticks, but up to a constant:
      axis.tickmode = "array";
      let selected = subsample(filteredSettings, maxTicks, firstHighlightedSample);

      axis.tickvals = selected.map((s) => s.name);
      axis.ticktext = selected.map((s) => "<span style='color:" + (s.highlight ?? "#ccc") + "'>" + s.name + "</span>");
    }
  }

  firstHighlightedSample(sampleSettings) {
    let index = 0;
    let highlighted = sampleSettings.filter((s) => s.highlight);
    if (highlighted.length > 0) index = sampleSettings.findIndex((s) => s.highlight);
    return index;
  }
}

function initPlot(dump) {
  if (dump["plot_type"] === "xy_line") return new LinePlot(dump);
  if (dump["plot_type"] === "bar_graph") return new BarPlot(dump);
  if (dump["plot_type"] === "box") return new BoxPlot(dump);
  if (dump["plot_type"] === "scatter") return new ScatterPlot(dump);
  if (dump["plot_type"] === "heatmap") return new HeatmapPlot(dump);
  if (dump["plot_type"] === "violin") return new ViolinPlot(dump);
  console.log("Did not recognise plot type: " + dump["plot_type"]);
  return null;
}

let loadingWarning;

$(function () {
  // Show loading warning
  loadingWarning = $(".mqc_loading_warning").show();
});

callAfterDecompressed.push(function (mqc_plotdata) {
  mqc_plots = Object.fromEntries(Object.values(mqc_plotdata).map((data) => [data.anchor, initPlot(data)]));

  let shouldLoad = $(".hc-plot.not_loaded:visible");

  // Show plots on page load: either render, or show the "Show Plot" button
  shouldLoad.each(function () {
    let anchor = $(this).attr("id");
    let plot = mqc_plots[anchor];
    setTimeout(function () {
      // Deferring each plot call prevents browser from locking up
      if (plot.deferRender) {
        $("#" + anchor)
          .removeClass("not_loaded")
          .html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
      } else {
        renderPlot(anchor);
      }
      if ($(".hc-plot.not_loaded:visible").length === 0)
        // All plots loaded successfully (rendered or deferred with "Show Plot"), so hiding the warning
        $(".mqc_loading_warning").hide();
    }, 50);
  });

  // All plots loaded successfully, so hiding the warning
  if (shouldLoad.length === 0) loadingWarning.hide();

  // Render a plot when clicked (heavy plots are not automatically rendered by default)
  $("body").on("click", ".render_plot", function (e) {
    let plotAnchor = $(this).parent().attr("id");
    renderPlot(plotAnchor);
  });

  // Button "Render all plots" clicked, so rendering everything, and hiding the parent button object
  $("#mqc-render-all-plots").click(function () {
    $(".hc-plot").each(function () {
      renderPlot($(this).attr("id"));
    });
    $(this).parent().hide();
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
    let plotAnchor = $(this).data("plot-anchor");

    // Toggling flags
    mqc_plots[plotAnchor].pActive = !$(this).hasClass("active");
    $(this).toggleClass("active");

    if (mqc_plots[plotAnchor].rendered) {
      renderPlot(plotAnchor); // re-render
    }
  });

  // A "Log" button above a plot is clicked
  $("button.interactive-switch-group.log10-switch").click(function (e) {
    e.preventDefault();
    let plotAnchor = $(this).data("plot-anchor");

    // Toggling flags
    mqc_plots[plotAnchor].lActive = !$(this).hasClass("active");
    $(this).toggleClass("active");

    if (mqc_plots[plotAnchor].rendered) {
      renderPlot(plotAnchor); // re-render
    }
  });

  // Switch data source
  $(".interactive-switch-group.dataset-switch-group button").click(function (e) {
    e.preventDefault();
    let el = $(this);
    if (el.hasClass("active")) return;
    el.siblings("button.active").removeClass("active");
    el.addClass("active");
    let plotAnchor = el.data("plot-anchor");
    let activeDatasetIdx = mqc_plots[plotAnchor].activeDatasetIdx;
    let newDatasetIdx = el.data("datasetIndex");
    mqc_plots[plotAnchor].activeDatasetIdx = newDatasetIdx;
    if (activeDatasetIdx === newDatasetIdx) return;

    if (mqc_plots[plotAnchor].rendered) {
      renderPlot(plotAnchor); // re-render
    }
  });

  // Make divs height-draggable
  // http://jsfiddle.net/Lkwb86c8/
  $(".hc-plot:not(.no-handle)").each(function () {
    let el = $(this);
    if (!el.parent().hasClass("hc-plot-wrapper")) {
      el.wrap('<div class="hc-plot-wrapper"></div>');
    }
    if (!el.siblings().hasClass("hc-plot-handle")) {
      el.after('<div class="hc-plot-handle"><span></span><span></span><span></span></div>');
    }
    el.css({ height: "auto", top: 0, bottom: "6px", position: "absolute" });
  });

  $(".hc-plot-handle").on("mousedown", function (e) {
    let wrapper = $(this).parent();
    let plotAnchor = wrapper.children(".hc-plot")[0].id;
    let startHeight = wrapper.height();
    let pY = e.pageY;

    let doc = $(document);
    doc.on("mouseup", function () {
      // Clear listeners now that we've let go
      doc.off("mousemove");
      doc.off("mouseup");
      // Fire off a custom jQuery event for other javascript chunks to tie into
      // Bind to the plot div, which should have a custom ID
      $(wrapper.parent().find(".hc-plot, .beeswarm-plot")).trigger("mqc_plotresize");
    });

    $(document).on("mousemove", function (me) {
      let newHeight = startHeight + (me.pageY - pY) + 2; // 2 px for the border or something
      wrapper.css("height", newHeight);
      if (mqc_plots[plotAnchor] !== undefined) mqc_plots[plotAnchor].resize(newHeight - 7); // 7 is the height of the handle overlapping the plot wrapper
    });
  });
});

// Highlighting, hiding and renaming samples. Takes a list of samples, returns
// a list of objects: {"name": "new_name", "highlight": "#cccccc", "hidden": false}
function applyToolboxSettings(samples, plotAnchor) {
  // init object with default values
  let objects = samples.map((name) => ({ name: name, highlight: null, hidden: false }));

  // Rename samples
  if (window.mqc_rename_f_texts.length > 0) {
    objects.map((obj) => {
      for (let patternIdx = 0; patternIdx < window.mqc_rename_f_texts.length; patternIdx++) {
        let pattern = window.mqc_rename_f_texts[patternIdx];
        let new_text = window.mqc_rename_t_texts[patternIdx];
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
    let groupDiv = $("#" + plotAnchor).closest(".mqc_hcplot_plotgroup");
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
      return objects;
    }
  }

  // Return the object indexed by sample names
  return objects;
}

// Call to render any plot
function renderPlot(plotAnchor) {
  let plot = mqc_plots[plotAnchor];
  if (plot === undefined) return false;
  if (plot.datasets.length === 0) return false;

  let container = $("#" + plotAnchor);

  // When the plot was already rendered, it's faster to call Plotly.react, using the same signature
  // https://plotly.com/javascript/plotlyjs-function-reference/#plotlyreact
  let func;
  if (!plot.rendered) {
    func = Plotly.newPlot;
    plot.rendered = true;
    container.removeClass("not_rendered").removeClass("not_loaded").parent().find(".render_plot").remove();
    if ($(".hc-plot.not_rendered").length === 0) $("#mqc-warning-many-samples").hide();
  } else {
    func = Plotly.react;
  }

  let dataset = plot.datasets[plot.activeDatasetIdx];
  updateObject(plot.layout, dataset.layout, false);

  // Apply pct/log toggle states
  plot.axisControlledBySwitches.map((axis) => {
    // Setting range explicitly just for the bar plot:
    plot.layout[axis].type = "linear";
    let min = plot.layout[axis]["autorangeoptions"]["minallowed"];
    let max = plot.layout[axis]["autorangeoptions"]["maxallowed"];
    if (plot.pActive) {
      updateObject(plot.layout[axis], plot.pctAxisUpdate, false);
      min = dataset["pct_range"][axis]["min"];
      max = dataset["pct_range"][axis]["max"];
    }
    if (plot.lActive) {
      plot.layout[axis].type = "log";
      // otherwise Plotly will interpret the range as log10:
      min = min && min > 0 ? Math.log10(min) : null;
      max = min && max > 0 ? Math.log10(max) : null;
    }
    plot.layout[axis]["autorangeoptions"]["minallowed"] = min;
    plot.layout[axis]["autorangeoptions"]["maxallowed"] = max;
  });

  let traces = plot.buildTraces();
  if (traces.length > 0 && traces[0].constructor === Array) traces = [].concat.apply([], traces); // if list of lists, flatten
  if (traces.length === 0) {
    // All series hidden. Hide the graph.
    container.hide();
    return;
  }

  container.show();
  func(plotAnchor, traces, plot.layout, {
    responsive: true,
    displaylogo: false,
    displayModeBar: true,
    toImageButtonOptions: { filename: plotAnchor },
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

function updateObject(target, source, nullOnly = false) {
  // Iterate through all keys in the source object
  for (const key in source) {
    // Check if the value is not null
    // Check if the value is an object and not an array
    if (Array.isArray(source[key])) {
      // Recursively update the array
      if (!nullOnly || target[key] === undefined || target[key] === null) {
        target[key] = [...source[key]];
      }
    } else if (typeof source[key] === "object") {
      // If the target doesn't have this key, or it's not an object, initialize it
      if (!target[key] || typeof target[key] !== "object") {
        target[key] = {};
      }
      // Recursively update the object
      updateObject(target[key], source[key], nullOnly);
    } else {
      if (!nullOnly || target[key] === undefined || target[key] === null) {
        // Directly update the value
        target[key] = source[key];
      }
    }
  }
}
