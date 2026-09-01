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

// Shared number formatter for tooltips and data labels across every echarts plot type
// (bar/line/scatter/box/heatmap/violin). Generalizes the ad hoc
// `Number.isInteger(v) ? v : parseFloat(v.toFixed(2))` pattern already used throughout the
// default (Plotly) template's plots/*.js for the same purpose, so hovering a point shows
// e.g. "0.01" instead of raw float noise like "0.010144431533554423". Values very close to
// zero fall back to a few significant digits instead of toFixed(2), which would otherwise
// collapse them to "0.00".
function formatNumber(value) {
  if (typeof value !== "number" || !Number.isFinite(value)) return value;
  if (Number.isInteger(value)) return value;
  if (value === 0 || Math.abs(value) >= 0.01) return parseFloat(value.toFixed(2));
  return Number(value.toPrecision(3));
}

// Make formatNumber globally available (single shared helper, reused by every plots/*.js
// tooltip/label formatter instead of a per-file copy).
window.formatNumber = formatNumber;

// SI-prefix abbreviation for value/log axis tick labels (POLISH.md #12), matching
// Plotly's default axis style: large numbers get k/M/G/T suffixes at ~3 significant
// figures with trailing zeros trimmed (450000000 -> "450M", 1500000 -> "1.5M",
// 12000 -> "12k"); small numbers pass through with the same rounding but no suffix.
// `suffix` is an optional caller-supplied unit appended verbatim (e.g. "x" for coverage),
// combining with the SI abbreviation rather than the two fighting over the same tick.
//
// GOLDEN CROSS-LANGUAGE CONTRACT: kept in lockstep with `_si_axis_formatter_body()` in
// `multiqc/plots/echarts/converter.py` (same duplication pattern as violin's `kde()`
// pair), since that Python copy is the only one ever executed (SSR path, via a `__FN__`
// sentinel); this JS copy is what runs for the live interactive render, applied in
// `buildCurrentOption()` below, which always overwrites the sentinel object before
// `setOption()` so it never reaches ECharts as-is (same pattern violin.js uses for its
// yAxis formatter).
function formatAxisNumber(v, suffix) {
  suffix = suffix || "";
  if (typeof v !== "number" || !Number.isFinite(v)) return String(v);
  const sign = v < 0 ? "-" : "";
  const abs = Math.abs(v);
  const units = [
    [1e12, "T"],
    [1e9, "G"],
    [1e6, "M"],
    [1e3, "k"],
  ];
  for (const [threshold, unitSuffix] of units) {
    if (abs >= threshold) {
      // ponytail: no rollover guard for a value that rounds up into the next unit
      // (e.g. 999.6k -> "1000k" not "1M"); axis ticks are always the round numbers
      // ECharts itself chooses, so this boundary never occurs in practice.
      return sign + Number((abs / threshold).toPrecision(3)) + unitSuffix + suffix;
    }
  }
  return sign + Number(abs.toPrecision(3)) + suffix;
}

// Make formatAxisNumber globally available (used by buildCurrentOption() below and by
// any per-type buildSeries()/applyOptionOverrides() that needs the same axis formatting).
window.formatAxisNumber = formatAxisNumber;

// d3-format specifier subset this understands (same subset as `_TICKFORMAT_RE` in
// multiqc/plots/echarts/converter.py::_tickformat_axis_formatter_body, kept in lockstep
// with that Python regex so a tooltip and its axis tick honour the same spec): an
// optional leading `,` (thousands grouping), an optional `.N` precision, then a single
// type letter (`f`/`%`/`s`, or bare).
const _HOVERFORMAT_RE = /^(,)?(?:\.(\d+))?([a-zA-Z%]?)$/;

// Formats `value` per `spec` (a Plotly/d3-format `hoverformat` string, e.g. ",.0f" /
// ".1%" / ".2s"), mirroring how the Plotly template's per-type tooltip formatters read
// `layout.xaxis.hoverformat`/`layout.yaxis.hoverformat` (see e.g.
// templates/plotly/src/js/plots/bar.js). Falls back to the existing formatNumber()
// heuristic when `spec` is falsy OR unsupported/malformed -- a bad or missing hoverformat
// must never break a tooltip, so this never throws.
function formatWithHoverformat(value, spec) {
  if (typeof value !== "number" || !Number.isFinite(value)) return value;
  if (!spec) return formatNumber(value);
  try {
    let m = _HOVERFORMAT_RE.exec(spec);
    if (!m) return formatNumber(value);
    let comma = m[1];
    let prec = m[2] !== undefined ? parseInt(m[2], 10) : null;
    let type = m[3];

    if (type === "s") return formatAxisNumber(value, "");
    if (type === "%") {
      let decimals = prec !== null ? prec : 0;
      return (value * 100).toFixed(decimals) + "%";
    }
    if (type === "f" || type === "") {
      let decimals = prec !== null ? prec : 6;
      return value.toLocaleString("en-US", {
        minimumFractionDigits: decimals,
        maximumFractionDigits: decimals,
        useGrouping: !!comma,
      });
    }
    return formatNumber(value);
  } catch (e) {
    return formatNumber(value);
  }
}

// Make formatWithHoverformat globally available, reused by every plots-echarts/*.js
// tooltip formatter that has an axis (or, for violin, a per-metric header) hoverformat
// to honour.
window.formatWithHoverformat = formatWithHoverformat;

// Take ~`num` samples from values evenly, always include `start`. If `roundBin` is true,
// the bins will be rounded to the nearest integer, so the ticks will be evenly
// distributed, but the total number of ticks may be less than `num`.
// Lifted to module scope (was nested inside Plot.recalculateTicks) so both
// recalculateTicks() and Plot.sampleAxisLabel() below can share one implementation
// (DRY -- ticket asked for exactly one shared subsample, not a per-caller copy).
function subsample(values, num, start = 0, roundBin = true) {
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

  // Reads the static markArea/markLine payload for the active dataset's
  // `x_bands`/`y_bands`/`x_lines`/`y_lines` (precomputed once in Python by
  // multiqc/plots/echarts/converter.py::bands_and_lines, stashed in the option
  // skeleton under `_mqc.bandsLines`) and, if present, wraps it as a trailing,
  // invisible `{type: "line"}` series. renderPlot() always replaces `option.series`
  // wholesale, so this has to be re-appended by every buildSeries() call rather than
  // set once; markArea/markLine attach fine to a "line" series even when the rest of
  // the chart's series are "bar"/"scatter", since they share the same grid/axes.
  // Called by line.js/bar.js/scatter.js's buildSeries() (DRY: one implementation).
  bandsLinesSeries() {
    let bandsLines = this.echarts.datasets[this.activeDatasetIdx].layout["_mqc"]?.bandsLines;
    if (!bandsLines) return [];
    return [
      {
        type: "line",
        name: "",
        data: [],
        silent: true,
        showSymbol: false,
        tooltip: { show: false },
        ...bandsLines,
      },
    ];
  }

  afterPlotCreated() {}

  // Hook for a plot type to mutate the fully-assembled option in place before it's set
  // on the chart (e.g. attaching a tooltip formatter function, which can never live in
  // the JSON-safe serialized skeleton). Called at the end of buildCurrentOption(),
  // after series (and any type-specific axis data) are in place.
  applyOptionOverrides(option) {}

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

  // ECharts equivalent of recalculateTicks() above: ECharts axisLabel ignores Plotly's
  // tickvals/ticktext, so highlight -> tick-label colouring here is instead a rich-text
  // `axisLabel.formatter` plus a `rich` style map keyed per sample. Returns undefined
  // when no highlight is active ANYWHERE (global flag, matching the dim/grey rule used
  // everywhere else, item A of the toolbox-highlight parity work) so the caller can leave
  // the axis's normal theme-colored axisLabel completely untouched. Shared by bar.js/
  // box.js (one sample axis, applied via Plot._axisData.axisLabel in
  // buildCurrentOption() below) and heatmap.js (two sample axes, applied directly since
  // heatmap needs both at once).
  sampleAxisLabel(filteredSettings, maxTicks) {
    if (window.mqc_highlight_f_texts.length === 0) return undefined;

    let firstHighlightedSample = this.firstHighlightedSample(filteredSettings);
    let selected = subsample(filteredSettings, maxTicks, firstHighlightedSample);
    let selectedNames = new Set(selected.map((s) => s.name));

    // Rich-text style keys must be simple identifiers, not raw sample names (which can
    // contain spaces/punctuation ECharts' rich-text mini-language can't parse as a style
    // key), so index them instead of keying by name directly.
    let rich = {};
    let keyByName = {};
    filteredSettings.forEach((s, i) => {
      let key = "s" + i;
      rich[key] = { color: s.highlight ?? "#ccc" };
      keyByName[s.name] = key;
    });

    let formatter = (value) => {
      if (!selectedNames.has(value)) return ""; // subsampled out, same as Plotly
      let key = keyByName[value];
      return key === undefined ? value : "{" + key + "|" + value + "}";
    };

    return { formatter, rich };
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
  if (dump["plot_type"] === "x/y line") return new window.EchartsLinePlot(dump);
  if (dump["plot_type"] === "scatter plot") return new window.EchartsScatterPlot(dump);
  if (dump["plot_type"] === "heatmap") return new window.EchartsHeatmapPlot(dump);
  if (dump["plot_type"] === "box plot") return new window.EchartsBoxPlot(dump);
  if (dump["plot_type"] === "violin plot") return new window.EchartsViolinPlot(dump);
  if (dump["plot_type"] === "seqcontent") return new window.EchartsSeqContentPlot(dump);
  // Every plot type is supported. An unknown one is a bug in the Python serializer
  // (multiqc/plots/echarts/__init__.py raises NotImplementedError for these), so fail
  // loudly here rather than rendering a placeholder card.
  throw new Error("Unknown ECharts plot type: " + dump["plot_type"]);
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

  // 1-5. Build the complete option (skeleton + theme + log/pct + toolbox-aware series).
  let option = buildCurrentOption(plot);
  if (!option) {
    // All series hidden. Hide the graph, same as plotting.js.
    container.hide();
    return;
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
  let isFirstChartInit = plot.chart === null;
  plot.chart ??= echarts.init(el, null, { renderer: plot.echarts.renderer });
  plot.chart.setOption(option, { notMerge: true });

  // Plotly-style click+drag box zoom (POLISH.md #17): activate ECharts' toolbox
  // dataZoom "brush select" cursor mode so the user can drag-to-zoom immediately,
  // without first clicking a toolbar button (there is none visible; see converter.py's
  // `zoom_option`, which pushes the toolbox icon row off-canvas). `notMerge: true`
  // above resets the chart's cursor mode on every render, so this has to be re-applied
  // every time, not just on first init. Plots with no zoomable axis (violin) never get
  // a `toolbox.feature.dataZoom` in their option, so this is a no-op for them.
  if (option.toolbox && option.toolbox.feature && option.toolbox.feature.dataZoom) {
    plot.chart.dispatchAction({
      type: "takeGlobalCursor",
      key: "dataZoomSelect",
      dataZoomSelectActive: true,
    });
  }

  if (isFirstChartInit) {
    // Double-click to reset zoom (Plotly behavior): reset every dataZoom component to
    // its full range. Harmless (no-op) on a chart with no dataZoom component (violin).
    // Registered once per chart instance, not per render, since ECharts event
    // listeners persist across setOption calls. Must be `chart.getZr().on(...)`
    // (zrender's own event bus), not `chart.on(...)` (ECharts' higher-level event bus,
    // which never fires "dblclick": confirmed empirically, it's not one of the named
    // events ECharts re-emits there, only zrender sees the raw DOM-level double-click).
    plot.chart.getZr().on("dblclick", () => {
      plot.chart.dispatchAction({ type: "dataZoom", start: 0, end: 100 });
    });
  }

  if (isFirstChartInit && plot.square) {
    // Unlike Plotly (which honors `layout.width` directly on `newPlot`, forcing a square
    // figure from the first render), ECharts always sizes its canvas to the container's
    // current CSS size (full width here). Force the square aspect once up front, reusing
    // the same width-from-height calc the drag-resize handle already applies via
    // `Plot.resize()`/`EchartsHeatmapPlot.resize()` (plotting-shared.js's mousedown
    // handler on `.hc-plot-handle`).
    plot.resize(plot.layout.height);
  }

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

// Build the complete ECharts option for the plot's currently active dataset: skeleton
// clone + theme colors + log/pct switch state + toolbox-aware series (steps 2-5 of
// renderPlot() above). Shared by renderPlot() and the export helper (toolbox-export.js)
// so the exported image can never drift from what's on screen. Returns null when every
// series is hidden (toolbox "hide all"), same signal renderPlot uses to hide the container.
function buildCurrentOption(plot) {
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

  // 3a. Theme the Plotly-style crosshair guide lines (POLISH.md #13): dashed lines in
  // the theme's spike color, label box matching the tooltip's own colors. The `type:
  // "cross"`/`show`/`triggerOn` structure itself is set once, for every plot type, in
  // the shared skeleton (converter.py's `convert_layout`); violin.py overrides it to
  // `{show: false}` there (misleading axis values), so the `option.axisPointer?.show`
  // guard below skips theming for violin, matching its disabled state.
  const axisPointerTheme = {
    crossStyle: { color: colors.spike_color, type: "dashed" },
    label: {
      backgroundColor: colors.hoverlabel_bgcolor,
      color: colors.hoverlabel_fontcolor,
      borderColor: colors.hoverlabel_bordercolor,
      borderWidth: 1,
    },
  };
  if (option.tooltip?.axisPointer) {
    option.tooltip.axisPointer = { ...option.tooltip.axisPointer, ...axisPointerTheme };
  }
  if (option.axisPointer?.show) {
    option.axisPointer = { ...option.axisPointer, ...axisPointerTheme };
  }

  // 3b. Carry `config.plot_font_family` (report config, read from mqc_config -- see
  // echarts template's head.html override) onto the option so ECharts uses it, matching
  // the SSR path's single-family constraint (theme.py FONT_FAMILY). Left unset (ECharts
  // default) when the user hasn't configured a custom font.
  if (window.mqc_config && window.mqc_config.plot_font_family) {
    option.textStyle = { ...option.textStyle, fontFamily: window.mqc_config.plot_font_family };
  }

  // 3c. SI-abbreviate value/log axis tick labels (POLISH.md #12), replacing the
  // JSON-safe `__FN__` sentinel the skeleton carries for these axes (converter.py's
  // `_convert_axis`) with a real live function; category axes (sample names) are left
  // alone. Runs before the log/pct switch below so a percent axis's "{value}%" formatter
  // (set unconditionally there) always wins over SI abbreviation.
  [
    ["xAxis", "xaxis"],
    ["yAxis", "yaxis"],
  ].forEach(([echartsAxisName, plotlyAxisName]) => {
    let axis = option[echartsAxisName];
    if (!axis || axis.type === "category") return;
    let suffix = plot.layout?.[plotlyAxisName]?.ticksuffix || "";
    axis.axisLabel = { ...axis.axisLabel, formatter: (v) => formatAxisNumber(v, suffix) };
  });

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
  if (!series || series.length === 0) return null;
  // Plotly-style crosshair mouse cursor (POLISH.md #13), not ECharts' default hand/
  // pointer. One place for every plot type, rather than a per-type buildSeries() edit:
  // ECharts' own default is `series.cursor: "pointer"`, so an explicit per-series value
  // already set (there is none, currently) would still win over this.
  series.forEach((s) => {
    if (s.cursor === undefined) s.cursor = "crosshair";
  });
  option.series = series;
  // Category-axis `data` is toolbox-dependent (sample names get hidden/renamed), so it's
  // never in the skeleton; each plot class stashes { axis: "xAxis"|"yAxis", data: [...] }
  // on itself during buildSeries() when (and only when) its category axis is populated
  // at runtime (bar: always yAxis; line: xAxis, only for `pconfig.categories` plots).
  if (plot._axisData && option[plot._axisData.axis]) {
    option[plot._axisData.axis].data = plot._axisData.data;
    // Highlight -> tick-label colouring (item B): sampleAxisLabel() returns undefined
    // when no highlight is active, so the axis's normal theme-colored axisLabel (set
    // above in step 3) is left alone; otherwise its formatter/rich win over the flat
    // theme color, since rich-text per-sample colors are more specific.
    if (plot._axisData.axisLabel) {
      option[plot._axisData.axis].axisLabel = {
        ...option[plot._axisData.axis].axisLabel,
        ...plot._axisData.axisLabel,
      };
    }
  }

  plot.applyOptionOverrides(option);

  return option;
}

// Make buildCurrentOption available globally (used by toolbox-export.js)
window.buildCurrentOption = buildCurrentOption;

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
