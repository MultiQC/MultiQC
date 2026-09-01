// Import Bootstrap 5 JavaScript
// Tree-shaking should remove unused JS
import * as bootstrap from "bootstrap";

// Make Bootstrap available globally for compatibility
window.bootstrap = bootstrap;

// Import MultiQC JavaScript modules in the correct order.
// This mirrors the Plotly template's import order, with the Plotly-specific plotting.js
// swapped for our echarts-plotting.js and every plot type routed through a plots-echarts/
// builder. Order matters: render.js must be last, and echarts-plotting.js (which defines
// window.initPlot/window.renderPlot and window.Plot) must run before any plots/*.js that
// extend window.Plot.

// Core functionality first (engine-neutral, reused from the default template)
import "./decompress.js";
import "./multiqc.js";
import "./flat.js";
import "./plotting-shared.js";

// ECharts plotting engine (replaces default's plotting.js)
import "./echarts-plotting.js";

// "tables.js" is a fork of the default template's file (Task 2.3): the general stats
// table's "plot a column as a scatter" modal calls Plotly.newPlot() directly (not through
// initPlot/renderPlot), so it needed its own ECharts rewrite. See the header comment in
// ./tables.js for the copy-drift risk; everything else in that file is unchanged.
import "./tables.js";
import "./doi.js";
import "./statuses.js";

// Color mode toggle functionality
import "./color-mode.js";

// Toolbox modules (engine-neutral, reused from the default template)
import "./toolbox/constants.js";
import "./toolbox/utils.js";
import "./toolbox/filters.js";
import "./toolbox/highlights.js";
import "./toolbox/rename.js";
import "./toolbox/hide.js";
// ECharts-specific image export (PNG/SVG): replaces default's Plotly.toImage-based
// toolbox/export.js with an ECharts equivalent (Task 0.7). Data export (csv/tsv/json)
// is unchanged, inherited via each plot's exportData().
import "./toolbox-export.js";
import "./toolbox/ai.js";
import "./toolbox/save-load.js";
import "./toolbox/citations.js";
import "./toolbox/help.js";
import "./toolbox.js";

// Plot types. All seven (bar, line, scatter, heatmap, box, violin, seqcontent) are ported
// and imported below; an unknown type raises loudly in window.initPlot rather than
// rendering a placeholder (see echarts-plotting.js).
// Default's plots/*.js define `class XPlot extends Plot` using the bare global `Plot`
// identifier: because echarts-plotting.js (imported above) runs first and does
// `window.Plot = Plot`, that bare reference resolves to OUR echarts Plot base class,
// not Plotly's. We import each default plot file first to get its
// prepData/exportData/formatDatasetForAiPrompt, then our echarts subclass extends
// window.BarPlot / window.LinePlot / window.ScatterPlot.
import "./plots/bar.js";
import "./plots-echarts/bar.js";
import "./plots/line.js";
import "./plots-echarts/line.js";
import "./plots/scatter.js";
import "./plots-echarts/scatter.js";
// "heatmap" (Phase 1, Task 1.3) is STANDALONE: it extends window.Plot directly and does
// NOT import the default template's plots/heatmap.js. That file has a top-level
// `$(function () {...})` handler (zmin/zmax range sliders + cluster toggle) that calls
// `Plotly.restyle`, which would run at load time and crash under this template (no
// `Plotly` global). `./plots/heatmap.js` copies the prepData()/exportData() field access
// it needs and re-implements both handlers for ECharts instead.
import "./plots-echarts/heatmap.js";
// "box" (Phase 1, Task 1.4): same extends-the-default-class strategy as bar/line/scatter
// above (default box.js has no top-level Plotly handler, only an engine-neutral
// sort-toggle click handler reused as-is).
import "./plots/box.js";
import "./plots-echarts/box.js";
// "violin" (Phase 2, Task 2.2) is STANDALONE (extends window.Plot directly, does NOT
// import the default template's plots/violin.js): that class's buildTraces() builds one
// Plotly subplot per metric row, which has no ECharts equivalent (a single value-axis
// grid is used instead). See ./plots/violin.js for the full rationale.
import "./plots-echarts/violin.js";
// "seqcontent" (Phase 2, Task 2.2): same extends-the-default-class strategy as
// bar/line/scatter/box above (default seqcontent.js has no top-level Plotly handler).
import "./plots/seqcontent.js";
import "./plots-echarts/seqcontent.js";

// AI features
import "./ai-helpers.js";
import "./ai.js";

// Render script should be last as it initializes everything
import "./render.js";
