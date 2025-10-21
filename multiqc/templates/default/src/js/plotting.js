////////////////////////////////////////////////
// Plotly Plotting Code
////////////////////////////////////////////////

// Global plot data variable. Accessed in many other JavaScript files.
window.mqc_plots = {};

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
    this.plotType = dump["plot_type"];
    this.pconfig = dump["pconfig"];
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

  afterPlotCreated() {}

  plotAiHeader(view) {
    let prompt = "Plot type: " + this.plotType + "\n";
    return prompt;
  }

  formatDatasetForAiPrompt(dataset) {
    return "";
  }

  formatForAiPrompt(view) {
    // Prepare data to be sent to the LLM. LLM doesn't need things like colors, etc.
    let result = this.plotAiHeader(view) + "\n\n";

    if (this.datasets.length === 1) {
      return result + this.formatDatasetForAiPrompt(this.datasets[0]);
    }

    for (let dataset of this.datasets) {
      let formattedDataset = this.formatDatasetForAiPrompt(dataset);
      if (!formattedDataset) continue;
      result += "### " + dataset.label + "\n";
      result += "\n";
      result += formattedDataset;
      result += "\n\n";
    }

    return result;
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

// Make Plot class available globally for plot classes to extend
window.Plot = Plot;

function initPlot(dump) {
  if (dump["plot_type"] === "bar plot") return new BarPlot(dump);
  if (dump["plot_type"] === "x/y line") return new LinePlot(dump);
  if (dump["plot_type"] === "box plot") return new BoxPlot(dump);
  if (dump["plot_type"] === "scatter plot") return new ScatterPlot(dump);
  if (dump["plot_type"] === "violin plot") return new ViolinPlot(dump);
  if (dump["plot_type"] === "heatmap") return new HeatmapPlot(dump);
  console.log("Did not recognise plot type: " + dump["plot_type"]);
  return null;
}

let loadingWarning;

$(function () {
  // Show loading warning
  loadingWarning = $(".mqc_loading_warning").show();
});

window.callAfterDecompressed.push(function (mqc_plotdata) {
  window.mqc_plots = Object.fromEntries(Object.values(mqc_plotdata).map((data) => [data.anchor, initPlot(data)]));

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
          .html('<button class="btn btn-outline-secondary btn-lg render_plot">Show plot</button>');
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

function getPseudonym(sampleName) {
  // If anonymization is enabled for AI prompts, use the aiPseudonymMap built in Python runtime,
  // which defines pseudonyms for original sample names (before renaming)

  // See if toolbox switch is on
  const anonymizationEnabled = getStoredSampleAnonymizationEnabled();
  if (!anonymizationEnabled) return undefined;

  // aiPseudonymMap is built in Python runtime and defined in head.html
  if (!aiPseudonymMap) return undefined;

  // Check the map. Exact match?
  if (aiPseudonymMap[sampleName]) return aiPseudonymMap[sampleName];

  // Try replacing partial matches for cases like sample="SAMPLE1-SAMPLE2"
  // Start with the longest original name to avoid situations when one sample is a prefix of another
  let sortedOriginals = Object.keys(aiPseudonymMap).sort((a, b) => b.length - a.length);
  let result = sampleName;
  for (let original of sortedOriginals) {
    if (result.includes(original)) {
      result = result.replace(original, aiPseudonymMap[original]);
    }
  }
  return result;
}

class Sample {
  constructor(name) {
    this.originalName = name;
    this.name = name;
    this.highlight = null;
    this.hidden = false;
    this.pseudonym = getPseudonym(name);
  }
}

// Highlighting, hiding and renaming samples. Takes a list of samples, returns
// a list of objects: {"name": "new_name", "highlight": "#cccccc", "hidden": false}
function applyToolboxSettings(samples, plotAnchor) {
  // init object with default values, apply pseudonymization
  let objects = samples.map((name) => new Sample(name));

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
      const alert = `
      <div class="samples-hidden-warning alert alert-warning">
        ⚠ <strong>Warning:</strong> ${nHidden} samples hidden.
        <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose('#mqc_hidesamples', true); return false;">See toolbox.</a>
      </div>`;
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

// Make functions available globally
window.applyToolboxSettings = applyToolboxSettings;

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

  // Apply theme-aware colors for hover elements and text
  const colors = getPlotlyThemeColors();

  // Main font color
  if (plot.layout.font) {
    plot.layout.font.color = colors.textcolor;
  }

  // Title color
  if (plot.layout.title) {
    plot.layout.title.font ??= {};
    plot.layout.title.font.color = colors.textcolor;
  }

  // Legend color
  if (plot.layout.legend) {
    plot.layout.legend.font ??= {};
    plot.layout.legend.font.color = colors.textcolor;
  }

  // Hover label
  plot.layout.hoverlabel = plot.layout.hoverlabel || {};
  plot.layout.hoverlabel.bgcolor = colors.hoverlabel_bgcolor;
  plot.layout.hoverlabel.bordercolor = colors.hoverlabel_bordercolor;
  plot.layout.hoverlabel.font = plot.layout.hoverlabel.font || {};
  plot.layout.hoverlabel.font.color = colors.hoverlabel_fontcolor;

  // Set spike line colors for hover crosshairs (only if spikes are explicitly enabled)
  if (plot.layout.xaxis && plot.layout.xaxis.showspikes === true) {
    plot.layout.xaxis.spikecolor = colors.spike_color;
  }
  if (plot.layout.yaxis && plot.layout.yaxis.showspikes === true) {
    plot.layout.yaxis.spikecolor = colors.spike_color;
  }

  // Apply theme colors to all axes (for plots with multiple axes like violin plots)
  for (let i = 1; i <= 20; i++) {
    let xaxis = "xaxis" + (i === 1 ? "" : i);
    let yaxis = "yaxis" + (i === 1 ? "" : i);

    if (plot.layout[xaxis]) {
      plot.layout[xaxis].gridcolor = colors.gridcolor;
      plot.layout[xaxis].zerolinecolor = colors.zerolinecolor;
      plot.layout[xaxis].color = colors.axiscolor;
      plot.layout[xaxis].tickfont ??= {};
      plot.layout[xaxis].tickfont.color = colors.tickcolor;
      if (plot.layout[xaxis].showspikes === true) {
        plot.layout[xaxis].spikecolor = colors.spike_color;
      }
    }

    if (plot.layout[yaxis]) {
      plot.layout[yaxis].gridcolor = colors.gridcolor;
      plot.layout[yaxis].zerolinecolor = colors.zerolinecolor;
      plot.layout[yaxis].color = colors.axiscolor;
      plot.layout[yaxis].tickfont ??= {};
      plot.layout[yaxis].tickfont.color = colors.tickcolor;
      if (plot.layout[yaxis].showspikes === true) {
        plot.layout[yaxis].spikecolor = colors.spike_color;
      }
    }
  }

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

  // Set colorbar colors for heatmaps
  traces.forEach((trace) => {
    if (trace.type === "heatmap" && trace.showscale) {
      trace.colorbar = trace.colorbar || {};
      trace.colorbar.tickfont = trace.colorbar.tickfont || {};
      trace.colorbar.tickfont.color = colors.tickcolor;
      trace.colorbar.titlefont = trace.colorbar.titlefont || {};
      trace.colorbar.titlefont.color = colors.textcolor;
    }
  });

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

// Make renderPlot available globally
window.renderPlot = renderPlot;

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

// Make updateObject globally available
window.updateObject = updateObject;

// Function to get theme-aware colors for Plotly graphs
function getPlotlyThemeColors() {
  const isDark = document.documentElement.getAttribute("data-bs-theme") === "dark";

  if (isDark) {
    return {
      paper_bgcolor: "rgba(0,0,0,0)", // transparent
      plot_bgcolor: "rgba(0,0,0,0)", // transparent
      gridcolor: "rgba(180,180,180,0.25)", // lighter gray for dark mode
      zerolinecolor: "rgba(180,180,180,0.3)",
      axiscolor: "rgba(200,200,200,1)", // lighter gray for axis labels
      tickcolor: "rgba(220,220,220,1)", // even lighter for tick text
      textcolor: "rgba(220,220,220,1)", // text color for titles and legends
      modebar_color: "rgba(200, 200, 200, 0.5)",
      modebar_activecolor: "rgba(220, 220, 220, 1)",
      hoverlabel_bgcolor: "rgba(40,40,40,1)", // dark background for tooltips (opaque)
      hoverlabel_bordercolor: "rgba(100,100,100,1)", // gray border for tooltips
      hoverlabel_fontcolor: "rgba(220,220,220,1)", // light text for dark mode
      spike_color: "rgba(220,220,220,1)", // light spike line color
    };
  } else {
    return {
      paper_bgcolor: "rgba(0,0,0,0)", // transparent
      plot_bgcolor: "rgba(0,0,0,0)", // transparent
      gridcolor: "rgba(128,128,128,0.15)", // darker gray for light mode
      zerolinecolor: "rgba(128,128,128,0.2)",
      axiscolor: "rgba(100,100,100,1)", // darker gray for axis labels
      tickcolor: "rgba(80,80,80,1)", // even darker for tick text
      textcolor: "rgba(60,60,60,1)", // text color for titles and legends
      modebar_color: "rgba(100, 100, 100, 0.5)",
      modebar_activecolor: "rgba(80, 80, 80, 1)",
      hoverlabel_bgcolor: "rgba(255,255,255,1)", // white background for tooltips (opaque)
      hoverlabel_bordercolor: "rgba(100,100,100,1)", // gray border
      hoverlabel_fontcolor: "rgba(30,30,30,1)", // dark text
      spike_color: "rgba(80,80,80,1)", // dark spike line color
    };
  }
}

// Make getPlotlyThemeColors globally available
window.getPlotlyThemeColors = getPlotlyThemeColors;

// Function to update all rendered Plotly graphs with new theme colors
function updatePlotlyTheme() {
  const colors = getPlotlyThemeColors();

  // Update all rendered plots
  $(".hc-plot:not(.not_rendered)").each(function () {
    const anchor = $(this).attr("id");
    const plot = mqc_plots[anchor];
    if (!plot || !plot.rendered) return;

    // Update layout colors
    const layoutUpdate = {
      paper_bgcolor: colors.paper_bgcolor,
      plot_bgcolor: colors.plot_bgcolor,
      "font.color": colors.textcolor, // main font color for titles and legends
      "title.font.color": colors.textcolor, // title color
      "legend.font.color": colors.textcolor, // legend text color
      "xaxis.gridcolor": colors.gridcolor,
      "xaxis.zerolinecolor": colors.zerolinecolor,
      "xaxis.color": colors.axiscolor,
      "xaxis.tickfont.color": colors.tickcolor,
      "xaxis.title.font.color": colors.textcolor, // axis title color
      "yaxis.gridcolor": colors.gridcolor,
      "yaxis.zerolinecolor": colors.zerolinecolor,
      "yaxis.color": colors.axiscolor,
      "yaxis.tickfont.color": colors.tickcolor,
      "yaxis.title.font.color": colors.textcolor, // axis title color
      "modebar.color": colors.modebar_color,
      "modebar.activecolor": colors.modebar_activecolor,
      "hoverlabel.bgcolor": colors.hoverlabel_bgcolor,
      "hoverlabel.bordercolor": colors.hoverlabel_bordercolor,
      "hoverlabel.font.color": colors.hoverlabel_fontcolor,
    };

    // Only update spike colors if spikes are enabled
    if (plot.layout.xaxis && plot.layout.xaxis.showspikes === true) {
      layoutUpdate["xaxis.spikecolor"] = colors.spike_color;
    }
    if (plot.layout.yaxis && plot.layout.yaxis.showspikes === true) {
      layoutUpdate["yaxis.spikecolor"] = colors.spike_color;
    }

    // For plots with multiple axes (like violin plots), update those too
    for (let i = 2; i <= 20; i++) {
      if (plot.layout["xaxis" + i]) {
        layoutUpdate["xaxis" + i + ".gridcolor"] = colors.gridcolor;
        layoutUpdate["xaxis" + i + ".zerolinecolor"] = colors.zerolinecolor;
        layoutUpdate["xaxis" + i + ".color"] = colors.axiscolor;
        layoutUpdate["xaxis" + i + ".tickfont.color"] = colors.tickcolor;
        layoutUpdate["xaxis" + i + ".title.font.color"] = colors.textcolor;
      }
      if (plot.layout["yaxis" + i]) {
        layoutUpdate["yaxis" + i + ".gridcolor"] = colors.gridcolor;
        layoutUpdate["yaxis" + i + ".zerolinecolor"] = colors.zerolinecolor;
        layoutUpdate["yaxis" + i + ".color"] = colors.axiscolor;
        layoutUpdate["yaxis" + i + ".tickfont.color"] = colors.tickcolor;
        layoutUpdate["yaxis" + i + ".title.font.color"] = colors.textcolor;
      }
    }

    // Apply the layout update
    Plotly.relayout(anchor, layoutUpdate);

    // Update heatmap colorbar colors and violin plot scatter marker colors
    const isDarkMode = document.documentElement.getAttribute("data-bs-theme") === "dark";
    const plotDiv = document.getElementById(anchor);
    if (plotDiv && plotDiv.data) {
      plotDiv.data.forEach((trace, idx) => {
        // Update heatmap colorbar colors
        if (trace.type === "heatmap" && trace.showscale) {
          Plotly.restyle(
            anchor,
            {
              "colorbar.tickfont.color": colors.tickcolor,
              "colorbar.titlefont.color": colors.textcolor,
            },
            [idx],
          );
        }
        // Update violin plot scatter marker colors
        if (trace.type === "scatter" && trace.marker && trace.marker.color) {
          const color = trace.marker.color;
          let newColor = null;
          if (isDarkMode) {
            // Light mode colors -> dark mode colors
            if (color === "#000000") newColor = "#ffffff";
            else if (color === "#0b79e6") newColor = "#5dade2";
          } else {
            // Dark mode colors -> light mode colors
            if (color === "#ffffff") newColor = "#000000";
            else if (color === "#5dade2") newColor = "#0b79e6";
          }
          if (newColor) {
            Plotly.restyle(anchor, { "marker.color": newColor }, [idx]);
          }

          // Update hoverlabel for scatter traces
          Plotly.restyle(
            anchor,
            {
              "hoverlabel.bgcolor": isDarkMode ? "rgba(40,40,40,1)" : "rgba(255,255,255,1)",
              "hoverlabel.font.color": isDarkMode ? "rgba(220,220,220,1)" : "rgba(30,30,30,1)",
            },
            [idx],
          );
        }
      });
    }
  });

  // Update table scatter plot if it exists and is rendered
  const tableScatterPlot = document.getElementById("table_scatter_plot");
  if (tableScatterPlot && !$(tableScatterPlot).hasClass("not_rendered")) {
    const layoutUpdate = {
      paper_bgcolor: colors.paper_bgcolor,
      plot_bgcolor: colors.plot_bgcolor,
      "font.color": colors.textcolor,
      "title.font.color": colors.textcolor,
      "xaxis.gridcolor": colors.gridcolor,
      "xaxis.zerolinecolor": colors.zerolinecolor,
      "xaxis.color": colors.axiscolor,
      "xaxis.tickfont.color": colors.tickcolor,
      "xaxis.title.font.color": colors.textcolor,
      "yaxis.gridcolor": colors.gridcolor,
      "yaxis.zerolinecolor": colors.zerolinecolor,
      "yaxis.color": colors.axiscolor,
      "yaxis.tickfont.color": colors.tickcolor,
      "yaxis.title.font.color": colors.textcolor,
      "hoverlabel.bgcolor": colors.hoverlabel_bgcolor,
      "hoverlabel.bordercolor": colors.hoverlabel_bordercolor,
      "hoverlabel.font.color": colors.hoverlabel_fontcolor,
    };
    Plotly.relayout("table_scatter_plot", layoutUpdate);
  }
}

// Set up theme change observer when DOM is ready
$(function () {
  const themeObserver = new MutationObserver((mutations) => {
    mutations.forEach((mutation) => {
      if (mutation.type === "attributes" && mutation.attributeName === "data-bs-theme") {
        updatePlotlyTheme();
      }
    });
  });

  // Observe the html element for data-bs-theme attribute changes
  themeObserver.observe(document.documentElement, {
    attributes: true,
    attributeFilter: ["data-bs-theme"],
  });
});
