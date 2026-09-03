// Tree-shaken ECharts bundle builder for MultiQC.
// Produces a minified IIFE bundle (global `echarts`) with the exact module
// list needed by the echarts plotting engine: bar/line/scatter/custom/
// heatmap/boxplot charts, tooltip/legend/grid/dataZoom/markLine/markArea/
// visualMap/title components, and both SVG + canvas renderers.
import { build } from "esbuild";
import { gzipSync } from "node:zlib";
import { readFileSync } from "node:fs";

const entry = "echarts-entry.js";
const outfile = "assets/js/packages/echarts-6.1.0.custom.min.js";

await build({
  entryPoints: [entry],
  bundle: true,
  minify: true,
  format: "iife",
  globalName: "echarts",
  outfile,
  legalComments: "none",
  target: ["es2018"],
});

const buf = readFileSync(outfile);
const gzip = gzipSync(buf, { level: 9 });
const kb = (n) => (n / 1024).toFixed(0) + " KB";
console.log(
  `${outfile}: ${buf.length.toLocaleString()} B (${kb(buf.length)}) raw, ${gzip.length.toLocaleString()} B (${kb(gzip.length)}) gzip`,
);
