// Import Bootstrap 5 JavaScript
// Tree-shaking should remove unused JS
import * as bootstrap from "bootstrap";

// Make Bootstrap available globally for compatibility
window.bootstrap = bootstrap;

// Import MultiQC JavaScript modules in the correct order
// Core functionality first
import "./decompress.js";
import "./multiqc.js";
import "./flat.js";
import "./plotting.js";
import "./tables.js";
import "./doi.js";
import "./statuses.js";

// Color mode toggle functionality
import "./color-mode.js";

// Toolbox modules
import "./toolbox/constants.js";
import "./toolbox/utils.js";
import "./toolbox/filters.js";
import "./toolbox/highlights.js";
import "./toolbox/rename.js";
import "./toolbox/hide.js";
import "./toolbox/export.js";
import "./toolbox/ai.js";
import "./toolbox/save-load.js";
import "./toolbox/citations.js";
import "./toolbox/help.js";
import "./toolbox.js";

// Plot types
import "./plots/bar.js";
import "./plots/box.js";
import "./plots/line.js";
import "./plots/scatter.js";
import "./plots/heatmap.js";
import "./plots/violin.js";

// AI features
import "./ai-helpers.js";
import "./ai.js";

// Render script should be last as it initializes everything
import "./render.js";
