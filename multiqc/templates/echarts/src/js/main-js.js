// Import Bootstrap 5 JavaScript
// Tree-shaking should remove unused JS
import * as bootstrap from "bootstrap";

// Make Bootstrap available globally for compatibility
window.bootstrap = bootstrap;

// Import MultiQC JavaScript modules in the correct order.
// This mirrors multiqc/templates/default/src/js/main-js.js exactly, with the
// Plotly-specific plotting.js swapped for our echarts-plotting.js, and (for
// Phase 0) only the bar plot type ported. Order matters: render.js must be
// last, and echarts-plotting.js (which defines window.initPlot/window.renderPlot
// and window.Plot) must run before any plots/*.js that extend window.Plot.

// Core functionality first (engine-neutral, reused from the default template)
import "../../../default/src/js/decompress.js";
import "../../../default/src/js/multiqc.js";
import "../../../default/src/js/flat.js";
import "../../../default/src/js/plotting-shared.js";

// ECharts plotting engine (replaces default's plotting.js)
import "./echarts-plotting.js";

import "../../../default/src/js/tables.js";
import "../../../default/src/js/doi.js";
import "../../../default/src/js/statuses.js";

// Color mode toggle functionality
import "../../../default/src/js/color-mode.js";

// Toolbox modules (engine-neutral, reused from the default template)
import "../../../default/src/js/toolbox/constants.js";
import "../../../default/src/js/toolbox/utils.js";
import "../../../default/src/js/toolbox/filters.js";
import "../../../default/src/js/toolbox/highlights.js";
import "../../../default/src/js/toolbox/rename.js";
import "../../../default/src/js/toolbox/hide.js";
// ECharts-specific image export (PNG/SVG): replaces default's Plotly.toImage-based
// toolbox/export.js with an ECharts equivalent (Task 0.7). Data export (csv/tsv/json)
// is unchanged, inherited via each plot's exportData().
import "./toolbox-export.js";
import "../../../default/src/js/toolbox/ai.js";
import "../../../default/src/js/toolbox/save-load.js";
import "../../../default/src/js/toolbox/citations.js";
import "../../../default/src/js/toolbox/help.js";
import "../../../default/src/js/toolbox.js";

// Plot types.
// "bar" (Phase 0) and "line" (Phase 1, Task 1.1) are ported; other types are not
// imported here (their window.initPlot entries throw "not supported yet" in
// echarts-plotting.js). Default's plots/*.js define `class XPlot extends Plot` using
// the bare global `Plot` identifier: because echarts-plotting.js (imported above) runs
// first and does `window.Plot = Plot`, that bare reference resolves to OUR echarts Plot
// base class, not Plotly's. We import each default plot file first to get its
// prepData/exportData/formatDatasetForAiPrompt, then our echarts subclass extends
// window.BarPlot / window.LinePlot.
import "../../../default/src/js/plots/bar.js";
import "./plots/bar.js";
import "../../../default/src/js/plots/line.js";
import "./plots/line.js";

// AI features
import "../../../default/src/js/ai-helpers.js";
import "../../../default/src/js/ai.js";

// Render script should be last as it initializes everything
import "../../../default/src/js/render.js";
