////////////////////////////////////////////////
// ECharts Plotting Code
//
// Analog of multiqc/templates/default/src/js/plotting.js, but for the Apache
// ECharts backend instead of Plotly. See "The ECharts model->JSON contract"
// in multiqc-echarts-exploration/BUILD_PLAN.md for the full contract this
// file implements (the 6-step renderPlot sequence, buildSeries()).
////////////////////////////////////////////////

// Function to get theme-aware colors for ECharts graphs.
// Same shape and values as getPlotlyThemeColors (plotting.js), keyed on data-bs-theme.
function getEchartsThemeColors() {
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
      hoverlabel_bgcolor: "rgba(255,255,255,1)", // white background for tooltips (opaque)
      hoverlabel_bordercolor: "rgba(100,100,100,1)", // gray border
      hoverlabel_fontcolor: "rgba(30,30,30,1)", // dark text
      spike_color: "rgba(80,80,80,1)", // dark spike line color
    };
  }
}

// Make getEchartsThemeColors globally available
window.getEchartsThemeColors = getEchartsThemeColors;

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
    // ECharts-specific: the serialized {renderer, datasets: [{layout: <option skeleton>}]}
    // payload (multiqc/plots/echarts/__init__.py::serialize), and the live chart instance
    // (created lazily on first render).
    this.echarts = dump["echarts"];
    this.chart = null;
  }

  activeDatasetSize() {
    throw new Error("activeDatasetSize() not implemented");
  }

  buildSeries() {
    throw new Error("buildSeries() not implemented");
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
    // Kept from the Plotly implementation verbatim: it only mutates the (Plotly-shaped)
    // layout object passed in, which ECharts ignores, so it's harmless to inherit.
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

  resize(newHeight, newWidth) {
    if (newHeight === null || newHeight === undefined) {
      console.error("Plot.resize: newHeight is " + newHeight);
      return;
    }
    this.layout.height = newHeight;

    if (newWidth !== null && newWidth !== undefined) {
      this.layout.width = newWidth;
    } else if (this.square) {
      // noinspection JSSuspiciousNameCombination
      this.layout.width = newHeight;
    } else {
      this.layout.width = null;
    }

    if (this.chart) {
      this.chart.resize({ height: newHeight, width: newWidth || "auto" });
    }
  }
}

// Make Plot class available globally for echarts plot classes to extend
window.Plot = Plot;

function initPlot(dump) {
  if (dump["plot_type"] === "bar plot") return new window.EchartsBarPlot(dump);
  // Not yet ported to ECharts (Phases 1-2 of the build plan). Rather than throwing (which
  // would crash the whole report render, since Plot.interactive_plot serializes every plot
  // on the page), build a lightweight placeholder plot; renderPlot() shows a visible,
  // actionable message for it instead of a chart.
  return new window.Plot(dump);
}

// Make initPlot available globally
window.initPlot = initPlot;

// Call to render any plot
function renderPlot(anchor) {
  let plot = mqc_plots[anchor];
  if (plot === undefined) return false;
  if (plot.datasets.length === 0) return false;

  let container = $("#" + anchor);
  let el = document.getElementById(anchor);

  if (plot.echarts && plot.echarts.unsupported) {
    el.innerHTML =
      '<div class="alert alert-secondary">The ECharts plotting engine does not yet support ' +
      plot.echarts.unsupported +
      " plots. Switch to the default template to view this plot.</div>";
    container.show();
    if (!plot.rendered) {
      plot.rendered = true;
      container.removeClass("not_rendered").removeClass("not_loaded").parent().find(".render_plot").remove();
      if ($(".hc-plot.not_rendered").length === 0) $("#mqc-warning-many-samples").hide();
    }
    return;
  }

  // 1. Pick the active dataset.
  let dataset = plot.datasets[plot.activeDatasetIdx];

  // 2. Clone the skeleton option (must stay JSON-safe: no functions live in it).
  let option = structuredClone(plot.echarts.datasets[plot.activeDatasetIdx].layout);

  // 3. Apply theme-aware colors.
  const colors = getEchartsThemeColors();
  option.backgroundColor = colors.paper_bgcolor;
  if (option.title) {
    option.title.textStyle = { ...option.title.textStyle, color: colors.textcolor };
  }
  if (option.legend) {
    option.legend.textStyle = { ...option.legend.textStyle, color: colors.textcolor };
  }
  ["xAxis", "yAxis"].forEach((axisKey) => {
    let axis = option[axisKey];
    if (!axis) return;
    axis.axisLine = { ...axis.axisLine, lineStyle: { ...axis.axisLine?.lineStyle, color: colors.axiscolor } };
    axis.axisLabel = { ...axis.axisLabel, color: colors.tickcolor };
    axis.nameTextStyle = { ...axis.nameTextStyle, color: colors.textcolor };
    axis.splitLine = { ...axis.splitLine, lineStyle: { ...axis.splitLine?.lineStyle, color: colors.gridcolor } };
  });
  if (option.tooltip) {
    option.tooltip.backgroundColor = colors.hoverlabel_bgcolor;
    option.tooltip.borderColor = colors.hoverlabel_bordercolor;
    option.tooltip.textStyle = { ...option.tooltip.textStyle, color: colors.hoverlabel_fontcolor };
  }

  // 4. Apply log/pct toggle states over the axes controlled by the plot's switches.
  plot.axisControlledBySwitches.forEach((axisName) => {
    let echartsAxisName = axisName === "xaxis" ? "xAxis" : "yAxis";
    let axis = option[echartsAxisName];
    if (!axis) return;
    if (plot.pActive) {
      let pctRange = (dataset["pct_range"] && dataset["pct_range"][axisName]) || {};
      axis.min = pctRange.min ?? 0;
      axis.max = pctRange.max ?? 100;
      axis.axisLabel = { ...axis.axisLabel, formatter: "{value}%" };
    }
    if (plot.lActive) {
      // ECharts takes raw min/max on log axes; unlike Plotly, no log10() conversion needed.
      axis.type = "log";
    }
  });

  // 5. Build series from the raw dataset fields (toolbox-aware).
  let series = plot.buildSeries();
  if (!series || series.length === 0) {
    // All series hidden. Hide the graph, same as plotting.js.
    container.hide();
    return;
  }
  option.series = series;
  if (plot._axisData && option.yAxis) {
    option.yAxis.data = plot._axisData;
  }

  container.show();

  // Class bookkeeping identical to the Plotly renderPlot's first-render handling.
  if (!plot.rendered) {
    plot.rendered = true;
    container.removeClass("not_rendered").removeClass("not_loaded").parent().find(".render_plot").remove();
    if ($(".hc-plot.not_rendered").length === 0) $("#mqc-warning-many-samples").hide();
  }

  // 6. Get-or-init the chart instance, then always setOption with notMerge so stale
  // series/axis state from a previous render never leaks through.
  plot.chart ??= echarts.init(el, null, { renderer: plot.echarts.renderer });
  plot.chart.setOption(option, { notMerge: true });

  if (!plot._resizeObserver) {
    plot._resizeObserver = new ResizeObserver(() => {
      if (plot.chart) plot.chart.resize();
    });
    plot._resizeObserver.observe(el);
  }

  plot.afterPlotCreated();
}

// Make renderPlot available globally
window.renderPlot = renderPlot;

// Re-render every already-rendered plot when the color theme changes, so the option is
// rebuilt from the skeleton with new theme colors (no echarts setTheme dependency).
$(function () {
  const themeObserver = new MutationObserver((mutations) => {
    mutations.forEach((mutation) => {
      if (mutation.type === "attributes" && mutation.attributeName === "data-bs-theme") {
        Object.keys(window.mqc_plots).forEach((anchor) => {
          if (window.mqc_plots[anchor].rendered) window.renderPlot(anchor);
        });
      }
    });
  });

  themeObserver.observe(document.documentElement, {
    attributes: true,
    attributeFilter: ["data-bs-theme"],
  });
});
