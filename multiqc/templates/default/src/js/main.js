// Import our custom CSS
import "../scss/main.scss";

// Import Bootstrap 5 JavaScript
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
