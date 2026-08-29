// ECharts violin plot: STANDALONE class extending window.Plot directly (Task 2.2).
// Unlike bar/line/scatter/box, this does NOT extend the default template's ViolinPlot:
// that class's buildTraces() builds one Plotly subplot (its own x/y axis pair) per
// metric row; this class does the ECharts equivalent itself, PLOTLY-STYLE PER-ROW
// SUBPLOTS: one ECharts `grid` (with its own `xAxis`/`yAxis` pair) per visible metric,
// stacked vertically, matching multiqc/plots/echarts/violin.py's own layout_option().
// Only the engine-neutral bits are ported here verbatim: `prepData()` field access
// (templates/default/src/js/plots/violin.js:18-43) and `exportData()` (:401-431).
//
// Rendering strategy (mirrors multiqc/plots/echarts/violin.py):
// - Each row is a `custom` series whose renderItem draws a closed KDE polygon in REAL
//   x-coordinates (each row has its own real-valued x-axis, POLISH.md #6) and a small
//   y-offset symmetric around 0 (the violin's density thickness; no row index baked in
//   any more, since each row is its own grid).
// - `plot.echarts.datasets[i].violins` ({metric: {poly, range}}) is precomputed by
//   Python from ALL samples. When no sample is hidden by the toolbox we reuse that
//   polygon directly verbatim (its y no longer depends on row position, so no
//   translation is needed even if a metric's row shifted); when a sample IS hidden we
//   recompute the KDE in JS from the visible values only (kde() below, a byte-for-byte
//   port of violin.py::kde -- see the golden fixture comment above it).
// - A beeswarm `scatter` series per row uses ECharts 6's native `jitter`.
// - `grid`/`xAxis`/`yAxis` are rebuilt wholesale every render (buildSeries(), applied in
//   applyOptionOverrides()) from the LIVE metric list, since it can differ from the
//   Python skeleton's SSR-time list (e.g. a table column toggled after page load); each
//   row gets a real, SI-formatted x-axis (POLISH.md #6) and a per-row `inside`/toolbox
//   dataZoom so click+drag zoom only affects that row (POLISH.md #8); the row's metric
//   title is drawn via the VALUE-AXIS TRICK y-axis (a single tick at y=0, the exact
//   vertical center of the row, POLISH.md #1).
// - Cross-row sample highlight: hovering a beeswarm point highlights that sample's
//   point on every other row via chart.on("mouseover")/dispatchAction("highlight").

// --- KDE: the cross-language contract (multiqc/plots/echarts/violin.py::kde) ---------
//
// GOLDEN VALUES (must match tests/test_plots_echarts.py::test_kde_golden_values and
// multiqc/plots/echarts/violin.py's own golden-value comment exactly):
//   values = [1.0, 2.0, 3.0, 4.0, 5.0]
//   xs     = [1.0, 3.0, 5.0]
//   kde(values, xs) == [0.15916497933387785, 0.1802710624663249, 0.15916497933387785]
//
// Epanechnikov-kernel density estimate with Silverman's rule-of-thumb bandwidth. Keep
// this formula exactly as written; do not "clean it up" without re-deriving the golden
// values in both languages.
function kde(values, xs) {
  const n = values.length;
  if (n === 0) throw new Error("kde: values must not be empty");
  const mean = values.reduce((a, b) => a + b, 0) / n;
  const variance = values.reduce((a, v) => a + (v - mean) * (v - mean), 0) / n;
  const sd = Math.sqrt(variance);
  const bandwidth = Math.max(1.06 * sd * Math.pow(n, -1 / 5), 1e-9);
  return xs.map((x) => {
    let density = 0;
    for (const v of values) {
      const u = (x - v) / bandwidth;
      if (Math.abs(u) <= 1) density += 0.75 * (1 - u * u);
    }
    return density / (n * bandwidth);
  });
}

// One-shot sanity check against the golden fixture above: logs loudly (not silently) if
// the port ever drifts from violin.py's formula. Cheap enough to run unconditionally at
// bundle load.
(function checkKdeGoldenValues() {
  const got = kde([1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 3.0, 5.0]);
  const want = [0.15916497933387785, 0.1802710624663249, 0.15916497933387785];
  const ok = got.length === want.length && got.every((v, i) => Math.abs(v - want[i]) < 1e-9);
  if (!ok) console.error("EchartsViolinPlot: JS kde() no longer matches violin.py golden values", got, want);
})();

const ROW_HEIGHT = 0.42;
const N_BINS = 60;
const RANGE_PAD = 0.15;

// Beeswarm point defaults (POLISH.md #9/FIX #9): mirrors multiqc/plots/violin.py's own
// Dataset.create scatter_trace_params (marker: {size: 4, color: ..., opacity: 1}, no
// marker line/border) and its default-template client-side dark-mode swap (black ->
// white, blue -> lighter blue when the viewer is in dark mode), so the echarts beeswarm
// reads the same as Plotly's in both themes. SCATTER_SIZE is the non-highlighted point
// diameter; SCATTER_HIGHLIGHT_SIZE (used only when toolbox highlighting is active)
// matches the default template's own highlighted-point size.
const SCATTER_SIZE = 4;
const SCATTER_HIGHLIGHT_SIZE = 8;
// Plotly's own choice: black when any metric has a custom color (so grey dots don't
// clash with a color-coded violin), else blue (so a plain grey violin's dots read as
// interactive); mirrors multiqc/plots/echarts/violin.py's _scatter_base_color.
const SCATTER_COLOR_COLORED = "#000000";
const SCATTER_COLOR_PLAIN = "#0b79e6";
// Dark-mode swap for the two base colors above (templates/default/src/js/plots/violin.js's
// own isDarkMode ? ... swap), applied only to the base (non-highlighted) color.
const SCATTER_COLOR_DARK_SWAP = { [SCATTER_COLOR_COLORED]: "#ffffff", [SCATTER_COLOR_PLAIN]: "#5dade2" };

// Inner Q1-Q3 box + median line drawn over each violin (POLISH.md #10b), sized as
// fractions of ROW_HEIGHT; mirrors multiqc/plots/echarts/violin.py's constants of the
// same name.
const BOX_HALF_HEIGHT = 0.12;
const MEDIAN_HALF_HEIGHT = 0.16;

// Fill opacity for the violin body (POLISH.md #10a); mirrors violin.py's FILL_ALPHA.
const FILL_ALPHA = 0.3;
// Box fill is more solid than the violin body so the Q1-Q3 range reads clearly against it.
const BOX_FILL_ALPHA = 0.55;

// --- Per-row grid geometry: the cross-language contract (mirrors
// multiqc/plots/echarts/violin.py's `_row_geometry`/ROW_PX/EXTRA_PX/BOTTOM_PX/
// ROW_GRID_FRACTION). Unlike the Python copy (whose actual SSR container height is
// outside its control), THIS copy also picks the container's CSS height directly (see
// buildSeries() below), so here the row proportions are exact, not just approximate:
// compact, Plotly-comparable rows (POLISH.md #4).
const ROW_PX = 42;
const EXTRA_PX = 55; // space reserved for the chart title/subtitle above the first row
const BOTTOM_PX = 10; // breathing room below the last row's x-axis tick labels
// Fraction of a row's slot handed to ECharts as the grid box (`containLabel: true` then
// shrinks the ACTUAL drawing area within that box to make room for the row's own x-axis
// tick labels, so this only needs to cover the tick labels plus a small gap to the next
// row, NOT also the violin height on top of that; mirrors violin.py's ROW_GRID_FRACTION).
const ROW_GRID_FRACTION = 0.86;

function rowGeometry(rowIdx, n) {
  const total = ROW_PX * n + EXTRA_PX + BOTTOM_PX;
  const topPct = (EXTRA_PX / total) * 100;
  const bottomPct = (BOTTOM_PX / total) * 100;
  const rowSlotPct = (100 - topPct - bottomPct) / n;
  const top = topPct + rowIdx * rowSlotPct;
  const height = rowSlotPct * ROW_GRID_FRACTION;
  return [`${top}%`, `${height}%`];
}

// Mirrors violin.py::_grid_left: shared left inset (px) for every row's grid, sized from
// the LONGEST visible title so every row's plot area starts at the same x, and computed
// (not left purely to `containLabel`'s own auto-measurement) so it's never smaller than
// this floor even in a context where text measurement is unreliable.
const TITLE_CHAR_PX = 6.6;
const MIN_GRID_LEFT = 60;
const MAX_GRID_LEFT = 280;

function gridLeft(titles) {
  if (titles.length === 0) return MIN_GRID_LEFT;
  const maxLen = arrMax(titles.map((t) => t.length));
  return Math.min(MAX_GRID_LEFT, Math.max(MIN_GRID_LEFT, maxLen * TITLE_CHAR_PX + 14));
}

// Linear-interpolation quantile (mirrors multiqc/plots/echarts/box.py::_quantile, also
// duplicated in this template's own box.js's quantile(); re-duplicated here rather than
// imported because each templates/echarts/src/js/plots/*.js file is its own ES module
// scope, same reasoning as this file's kde() mirror above).
function quantile(sortedValues, q) {
  const n = sortedValues.length;
  if (n === 1) return sortedValues[0];
  const h = (n - 1) * q;
  const lo = Math.floor(h);
  const hi = Math.min(lo + 1, n - 1);
  const frac = h - lo;
  return sortedValues[lo] + frac * (sortedValues[hi] - sortedValues[lo]);
}

// Reduce-based min/max: `Math.min(...arr)`/`Math.max(...arr)` blow the call stack on
// large sample counts (general stats tables can have thousands of rows).
function arrMin(arr) {
  return arr.reduce((a, b) => (b < a ? b : a), arr[0]);
}
function arrMax(arr) {
  return arr.reduce((a, b) => (b > a ? b : a), arr[0]);
}

// Mirrors violin.py::_metric_range: prefers the column's own configured axis range
// (header.xaxis.range, the exact [dmin, dmax] Plotly's per-row x-axis uses) so a row's
// x-axis matches Plotly's exactly (POLISH.md #10c/#6), falling back to a padded
// observed-value range when the column has no configured range.
function metricRange(header, values) {
  const range = header && header.xaxis && header.xaxis.range;
  if (range && range[1] > range[0]) return [range[0], range[1]];

  let lo = arrMin(values);
  let hi = arrMax(values);
  const span = hi - lo;
  lo -= span * RANGE_PAD;
  hi += span * RANGE_PAD;
  if (hi <= lo) {
    // All values identical (span == 0): center a small synthetic window on the value
    // instead of pinning it to the row's left edge.
    lo = values[0] - 0.5;
    hi = values[0] + 0.5;
  }
  return [lo, hi];
}

// Mirrors violin.py::_violin_polygon: closed KDE polygon in REAL x-space (each row has
// its own real-valued x-axis) and a small y-space symmetric around 0 (no row offset any
// more, since each row is its own grid).
function violinPolygon(values, lo, hi) {
  const span = hi - lo;
  if (arrMax(values) === arrMin(values)) {
    // Degenerate (zero-variance) metric: kde()'s Silverman bandwidth collapses to its
    // 1e-9 floor here, which would otherwise blow up into a single needle-thin,
    // full-height spike (POLISH.md #10d). Draw a flat row instead; the median tick
    // (see buildSeries()) still marks the value.
    const x = values[0];
    return [
      [x, 0],
      [x, 0],
    ];
  }
  const xs = Array.from({ length: N_BINS }, (_, i) => lo + (span * i) / (N_BINS - 1));
  const ys = kde(values, xs);
  const ymax = arrMax(ys) || 1.0;
  const top = xs.map((x, i) => [x, (ys[i] / ymax) * ROW_HEIGHT]);
  const bottom = xs.map((x, i) => [x, -(ys[i] / ymax) * ROW_HEIGHT]).reverse();
  return top.concat(bottom);
}

// Mirrors converter.py's _clean_title_segment (POLISH.md #7): titles/namespaces are
// authored with HTML entities (e.g. "Mapped &amp; paired") for Plotly's HTML-interpreting
// axis labels, but ECharts renders axis labels as plain text. DOMParser decodes entities
// and strips any stray tags in one pass (no HTML re-injection risk: textContent never
// executes markup, it only reads the parsed text).
function decodeAndStripHtml(text) {
  if (!text) return text;
  const doc = new DOMParser().parseFromString(text, "text/html");
  return (doc.body.textContent || "").replace(/\s+/g, " ").trim();
}

// Mirrors violin.py::_metric_title.
function metricTitle(header) {
  const title = header.namespace ? `${header.namespace}: ${header.title}` : header.title;
  return decodeAndStripHtml(title);
}

// Copied from templates/default/src/js/plots/violin.js's buildTraces() (normalizeColorToRGB
// inline helper): accepts either "r,g,b" or "#rrggbb", returns "r,g,b" or null.
function normalizeColorToRGB(color) {
  if (!color) return null;
  if (/^\d+,\s*\d+,\s*\d+$/.test(color)) return color;
  if (color.startsWith("#")) {
    const hex = color.replace("#", "");
    if (hex.length === 6) {
      const r = parseInt(hex.substr(0, 2), 16);
      const g = parseInt(hex.substr(2, 2), 16);
      const b = parseInt(hex.substr(4, 2), 16);
      return `${r},${g},${b}`;
    }
  }
  return null;
}

// "rgb(r,g,b)" -> "rgba(r,g,b,alpha)"; mirrors violin.py::_rgba.
function toRgba(rgb, alpha) {
  return rgb.replace("rgb(", "rgba(").replace(")", `,${alpha})`);
}

// Default violin body color (POLISH.md #5): matches Plotly's own general-stats violin
// default, fillcolor="#999999" (multiqc/plots/violin.py's Dataset.create); mirrors
// violin.py's _DEFAULT_STROKE/_DEFAULT_FILL/_DEFAULT_BOX_FILL. Only used when a metric
// has no configured header.color; a real per-metric color always takes priority below.
const DEFAULT_STROKE = "rgb(153,153,153)";

// Mirrors violin.py::_violin_colors: [fill, stroke, boxFill], falling back to the same
// defaults (kept as rgb(...), not hex, so toRgba() above applies uniformly).
function violinColors(header) {
  const rgb = header.color ? normalizeColorToRGB(header.color) : null;
  if (!rgb) return [toRgba(DEFAULT_STROKE, FILL_ALPHA), DEFAULT_STROKE, toRgba(DEFAULT_STROKE, BOX_FILL_ALPHA)];
  const stroke = `rgb(${rgb})`;
  return [toRgba(stroke, FILL_ALPHA), stroke, toRgba(stroke, BOX_FILL_ALPHA)];
}

// Min/max row label: delegates to the shared window.formatNumber (echarts-plotting.js)
// so every echarts plot type rounds numbers the same way.
function formatG(num) {
  return String(window.formatNumber(num));
}

// Draws the KDE polygon plus an inner Q1-Q3 box and median line (POLISH.md #10b),
// matching Plotly's box_visible/meanline violin config. q1/median/q3 are real values
// (same x-space as poly). Mirrors violin.py::_violin_render_series.
function makeViolinRenderItem(poly, q1, median, q3, fill, stroke, boxFill) {
  return function (params, api) {
    const pts = poly.map((p) => api.coord(p));
    const bTL = api.coord([q1, -BOX_HALF_HEIGHT]);
    const bBR = api.coord([q3, BOX_HALF_HEIGHT]);
    const mTop = api.coord([median, -MEDIAN_HALF_HEIGHT]);
    const mBot = api.coord([median, MEDIAN_HALF_HEIGHT]);
    return {
      type: "group",
      children: [
        { type: "polygon", shape: { points: pts }, style: { fill, stroke, lineWidth: 1 } },
        {
          type: "rect",
          shape: {
            x: Math.min(bTL[0], bBR[0]),
            y: Math.min(bTL[1], bBR[1]),
            width: Math.abs(bBR[0] - bTL[0]),
            height: Math.abs(bBR[1] - bTL[1]),
          },
          style: { fill: boxFill, stroke, lineWidth: 1 },
        },
        {
          type: "line",
          shape: { x1: mTop[0], y1: mTop[1], x2: mBot[0], y2: mBot[1] },
          style: { stroke, lineWidth: 2.5 },
        },
      ],
    };
  };
}

// loObs/hiObs are the violin's real observed extent; drawn at the row's vertical center
// (y=0). Mirrors violin.py::_annotation_render_series.
function makeAnnotationRenderItem(loObs, hiObs) {
  const loLabel = formatG(loObs);
  const hiLabel = formatG(hiObs);
  return function (params, api) {
    const L = api.coord([loObs, 0]);
    const R = api.coord([hiObs, 0]);
    return {
      type: "group",
      children: [
        {
          type: "text",
          style: {
            text: loLabel,
            x: L[0] - 6,
            y: L[1],
            textAlign: "right",
            textVerticalAlign: "middle",
            fontSize: 10,
            fill: "#888",
          },
        },
        {
          type: "text",
          style: {
            text: hiLabel,
            x: R[0] + 6,
            y: R[1],
            textAlign: "left",
            textVerticalAlign: "middle",
            fontSize: 10,
            fill: "#888",
          },
        },
      ],
    };
  };
}

class EchartsViolinPlot extends window.Plot {
  constructor(dump) {
    super(dump);
    this.tableAnchor = dump["table_anchor"];
    this.isDownsampled = dump["is_downsampled"];
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0;
    return this.datasets[this.activeDatasetIdx]["all_samples"].length;
  }

  // Ported field-for-field from templates/default/src/js/plots/violin.js:16-83 (same
  // toolbox hide/rename/highlight handling, same table-column-visibility check).
  prepData(dataset, keepHidden = false) {
    dataset = dataset ?? this.datasets[this.activeDatasetIdx];
    let metrics = dataset["metrics"];
    let headerByMetric = dataset["header_by_metric"];

    metrics = metrics.filter((metric) => {
      let header = headerByMetric[metric];
      let tableCheckbox = $(`#${this.tableAnchor}_config_modal_table .mqc_table_col_visible[value="${metric}"]`);
      let hiddenInTable = tableCheckbox.length > 0 && !tableCheckbox.is(":checked");
      return (header["hidden"] !== true || keepHidden) && !hiddenInTable;
    });

    let violinValuesBySampleByMetric = dataset["violin_value_by_sample_by_metric"];
    let scatterValuesBySampleByMetric = {};
    metrics.forEach((metric) => {
      let header = headerByMetric[metric];
      let scatterValuesBySample = {};
      if (header["show_points"]) {
        if (header["show_only_outliers"]) scatterValuesBySample = dataset["scatter_value_by_sample_by_metric"][metric];
        else scatterValuesBySample = violinValuesBySampleByMetric[metric];
      }
      scatterValuesBySampleByMetric[metric] = scatterValuesBySample;
    });

    let allSamples = dataset["all_samples"];
    let sampleSettings = applyToolboxSettings(allSamples);

    let someHidden = sampleSettings.filter((s) => s.hidden).length > 0;
    if (someHidden) {
      let filteredViolinValuesBySampleByMetric = {};
      let filteredScatterValuesBySampleByMetric = {};
      metrics.map((metric) => {
        filteredViolinValuesBySampleByMetric[metric] = {};
        Object.keys(violinValuesBySampleByMetric[metric]).map((sample) => {
          if (!sampleSettings[allSamples.indexOf(sample)].hidden)
            filteredViolinValuesBySampleByMetric[metric][sample] = violinValuesBySampleByMetric[metric][sample];
        });
      });
      violinValuesBySampleByMetric = filteredViolinValuesBySampleByMetric;
      // NOTE: templates/default/src/js/plots/violin.js has the same someHidden block,
      // but (a) when show_points && show_only_outliers it reads from
      // filteredScatterValuesBySampleByMetric[metric] before that metric's entry has
      // been populated (Object.keys(undefined) throws), and (b) it never reassigns the
      // filtered result back to the returned scatterValuesBySampleByMetric, so hidden
      // samples never actually get filtered out of the beeswarm there. Both are latent
      // bugs in the default template; harmless to fix here since this is our own copy,
      // not an edit to the default template. The values to filter are already correct
      // per metric (scatterValuesBySampleByMetric was built above from show_points/
      // show_only_outliers), so this just needs to drop hidden samples from it.
      metrics.forEach((metric) => {
        let scatterValuesBySample = scatterValuesBySampleByMetric[metric] || {};
        filteredScatterValuesBySampleByMetric[metric] = {};
        Object.keys(scatterValuesBySample).forEach((sample) => {
          if (!sampleSettings[allSamples.indexOf(sample)].hidden)
            filteredScatterValuesBySampleByMetric[metric][sample] = scatterValuesBySample[sample];
        });
      });
      scatterValuesBySampleByMetric = filteredScatterValuesBySampleByMetric;
    }

    return [
      metrics,
      headerByMetric,
      allSamples,
      sampleSettings,
      violinValuesBySampleByMetric,
      scatterValuesBySampleByMetric,
    ];
  }

  // Copied verbatim (same field access) from templates/default/src/js/plots/violin.js:401-431.
  exportData(format) {
    let [metrics, headerByMetric, allSamples, sampleSettings, violinValuesBySampleByMetric] = this.prepData();

    let sep = format === "tsv" ? "\t" : ",";
    let titles = metrics.map((metric) => headerByMetric[metric].title);
    titles = titles.map((title) => (title.includes(sep) ? `"${title}"` : title));
    let csv = "Sample" + sep + titles.join(sep) + "\n";
    for (let i = 0; i < allSamples.length; i++) {
      let sample = allSamples[i];
      if (sampleSettings[i].hidden) continue;
      csv += sample + sep;
      csv += metrics
        .map((metric) => {
          let val = violinValuesBySampleByMetric[metric][sample];
          if (val === undefined) val = ".";
          return val;
        })
        .join(sep);
      csv += "\n";
    }
    return csv;
  }

  // Builds: n_metric violin (KDE polygon + inner box + median) `custom` series, then
  // n_metric min/max-annotation `custom` series, then n_metric beeswarm `scatter`
  // series -- same ordering as the SSR path (multiqc/plots/echarts/violin.py::series()),
  // which is what the sample -> [seriesIndex, dataIndex] cross-highlight map below relies
  // on. Also (re)builds `this._grids`/`this._xAxis`/`this._yAxis`/`this._toolbox`/
  // `this._dataZoom` from the LIVE metric list (PLOTLY-STYLE PER-ROW SUBPLOTS), applied
  // onto the option in applyOptionOverrides() below.
  buildSeries() {
    let [
      metrics,
      headerByMetric,
      allSamples,
      sampleSettings,
      violinValuesBySampleByMetric,
      scatterValuesBySampleByMetric,
    ] = this.prepData();

    const someHidden = sampleSettings.some((s) => s.hidden);
    const highlightingEnabled = sampleSettings.some((s) => s.highlight !== null);

    const echartsDs = this.echarts.datasets[this.activeDatasetIdx];
    const violinsFromPython = echartsDs.violins || {};

    // Pre-filter to metrics that actually have numeric data left to draw (mirrors
    // violin.py's _visible_metrics: filtered ONCE so the row index / row count used for
    // series AND for grid/axis geometry never disagree, even though a metric can drop
    // out here for reasons Python's SSR-time filter didn't know about, e.g. every one of
    // its samples got hidden via the toolbox).
    const rows = [];
    metrics.forEach((metric) => {
      const header = headerByMetric[metric];
      const valuesBySample = violinValuesBySampleByMetric[metric] || {};
      const numericValues = Object.values(valuesBySample).filter((v) => typeof v === "number" && Number.isFinite(v));
      if (numericValues.length === 0) return;
      rows.push({ metric, header, numericValues });
    });
    const n = rows.length;

    const height = ROW_PX * Math.max(n, 1) + EXTRA_PX + BOTTOM_PX;
    const el = document.getElementById(this.anchor);
    if (el) el.style.height = height + "px";
    $("#" + this.anchor + "-wrapper").css("height", height + "px");

    if (n === 0) {
      this._grids = [];
      this._xAxis = [];
      this._yAxis = [];
      this._rowCount = 0;
      return [];
    }

    const colors = window.getEchartsThemeColors();
    const rowLeft = gridLeft(rows.map(({ header }) => metricTitle(header)));
    const violinSeries = [];
    const annotationSeries = [];
    const scatterSeries = [];
    const grids = [];
    const xAxis = [];
    const yAxis = [];

    // Beeswarm base point color (FIX #9): one choice for the whole dataset, mirroring
    // multiqc/plots/violin.py's Dataset.create (black when any metric has a custom
    // color, else blue) plus the default template's own dark-mode swap, so the two
    // engines' beeswarms look the same in either viewer theme.
    const isDarkMode = document.documentElement.getAttribute("data-bs-theme") === "dark";
    const anyColored = Object.values(headerByMetric).some((h) => h.color);
    let scatterBaseColor = anyColored ? SCATTER_COLOR_COLORED : SCATTER_COLOR_PLAIN;
    if (isDarkMode) scatterBaseColor = SCATTER_COLOR_DARK_SWAP[scatterBaseColor] || scatterBaseColor;

    // sample (display name) -> [{title, value, suffix}, ...] across every row (FIX #3):
    // built alongside each row's scatterData below, so the combined cross-row tooltip
    // (wired via option.tooltip.formatter in applyOptionOverrides()) can show one
    // sample's value in every metric row at once without re-walking the series data.
    const sampleTooltipData = {};

    rows.forEach(({ metric, header, numericValues }, rowIdx) => {
      let poly, lo, hi;
      const pre = violinsFromPython[metric];
      if (!someHidden && pre) {
        // Fast path: reuse the precomputed polygon verbatim. Its y no longer encodes a
        // row offset (each row is its own grid), so unlike the old design no translation
        // is needed even if this metric's row position shifted relative to Python's
        // serialization order (e.g. a metric hidden via the table column config).
        [lo, hi] = pre.range;
        poly = pre.poly;
      } else {
        // Slow path: samples are hidden, so the polygon must be recomputed from the
        // currently-visible values only (see the perf note in BUILD_PLAN.md Phase 2 risks).
        [lo, hi] = metricRange(header, numericValues);
        poly = violinPolygon(numericValues, lo, hi);
      }

      // Q1/median/Q3 are cheap (a sort + linear interpolation) compared to the KDE, so
      // unlike the polygon they're just always recomputed fresh here, fast path or not
      // (see the module docstring's point 2 note on this).
      const sorted = [...numericValues].sort((a, b) => a - b);
      const q1 = quantile(sorted, 0.25);
      const median = quantile(sorted, 0.5);
      const q3 = quantile(sorted, 0.75);
      const obsMin = arrMin(numericValues);
      const obsMax = arrMax(numericValues);

      const [fill, stroke, boxFill] = violinColors(header);
      violinSeries.push({
        type: "custom",
        coordinateSystem: "cartesian2d",
        xAxisIndex: rowIdx,
        yAxisIndex: rowIdx,
        data: [0],
        renderItem: makeViolinRenderItem(poly, q1, median, q3, fill, stroke, boxFill),
        silent: true,
        z: 1,
      });

      annotationSeries.push({
        type: "custom",
        coordinateSystem: "cartesian2d",
        xAxisIndex: rowIdx,
        yAxisIndex: rowIdx,
        data: [0],
        renderItem: makeAnnotationRenderItem(obsMin, obsMax),
        silent: true,
        z: 3,
      });

      const title = metricTitle(header);
      const suffix = (header.xaxis && header.xaxis.ticksuffix) || "";

      const scatterData = [];
      if (header.show_points) {
        const svBySample = scatterValuesBySampleByMetric[metric] || {};
        Object.entries(svBySample).forEach(([sample, value]) => {
          if (typeof value !== "number" || !Number.isFinite(value)) return;
          const state = sampleSettings[allSamples.indexOf(sample)];
          let color = scatterBaseColor;
          let size = SCATTER_SIZE;
          if (highlightingEnabled) {
            color = state?.highlight ?? "#ddd";
            size = state?.highlight != null ? SCATTER_HIGHLIGHT_SIZE : SCATTER_SIZE;
          }
          const name = state?.name ?? sample;
          scatterData.push({
            // Real value: no denormalization needed any more, this row's x-axis is
            // already real-valued (POLISH.md #6).
            value: [value, 0],
            name,
            // FIX #9: no marker border/opacity fade, matching Plotly's own solid,
            // borderless beeswarm dots (multiqc/plots/violin.py's scatter_trace_params).
            itemStyle: { color, opacity: 1 },
            symbolSize: size,
          });
          // FIX #3: this row's contribution to the sample's combined cross-row tooltip.
          if (!sampleTooltipData[name]) sampleTooltipData[name] = [];
          sampleTooltipData[name].push({ title, value, suffix });
        });
      }
      scatterSeries.push({
        type: "scatter",
        name: String(metric),
        xAxisIndex: rowIdx,
        yAxisIndex: rowIdx,
        data: scatterData,
        symbolSize: SCATTER_SIZE,
        jitter: 22,
        jitterOverlap: false,
        // FIX #2: no persistent per-point text; the sample's info only ever shows in the
        // hover tooltip (option.tooltip.formatter, applyOptionOverrides() below).
        label: { show: false },
        z: 2,
      });

      // Grid/axis geometry for this row (PLOTLY-STYLE PER-ROW SUBPLOTS): a real,
      // SI-formatted x-axis (POLISH.md #6), and a hidden value y-axis whose single tick
      // (at y=0, this row's exact vertical center) carries the metric title
      // (POLISH.md #1, the VALUE-AXIS TRICK, see multiqc/plots/echarts/violin.py's
      // module docstring), themed to match every other plot type (mirrors the
      // buildCurrentOption theme step in echarts-plotting.js, which never runs on these
      // axes since they're built fresh here, after that step already ran).
      const [top, gridHeight] = rowGeometry(rowIdx, n);
      grids.push({ top, height: gridHeight, left: rowLeft, right: 16, containLabel: true });
      xAxis.push({
        type: "value",
        gridIndex: rowIdx,
        min: lo,
        max: hi,
        axisLabel: { fontSize: 10, color: colors.tickcolor, formatter: (v) => window.formatAxisNumber(v, suffix) },
        axisLine: { show: true, lineStyle: { color: colors.axiscolor } },
        axisTick: { show: true, lineStyle: { color: colors.axiscolor } },
        splitLine: { show: false },
      });
      yAxis.push({
        type: "value",
        gridIndex: rowIdx,
        min: -0.5,
        max: 0.5,
        interval: 0.5,
        axisLabel: {
          fontSize: 12,
          align: "right",
          verticalAlign: "middle",
          color: colors.tickcolor,
          formatter: (v) => (Math.abs(v) < 1e-6 ? title : ""),
        },
        axisTick: { show: false },
        axisLine: { show: false },
        splitLine: { show: false },
      });
    });

    this._grids = grids;
    this._xAxis = xAxis;
    this._yAxis = yAxis;
    this._rowCount = n;

    // sample -> [[seriesIndex, dataIndex], ...] across every scatter series, for the
    // cross-row highlight wired in afterPlotCreated().
    const sampleIndexMap = {};
    const scatterOffset = violinSeries.length + annotationSeries.length;
    scatterSeries.forEach((series, si) => {
      series.data.forEach((item, di) => {
        if (!sampleIndexMap[item.name]) sampleIndexMap[item.name] = [];
        sampleIndexMap[item.name].push([scatterOffset + si, di]);
      });
    });
    this._sampleIndexMap = sampleIndexMap;
    this._sampleTooltipData = sampleTooltipData;

    return violinSeries.concat(annotationSeries, scatterSeries);
  }

  // grid/xAxis/yAxis are rebuilt wholesale from the live metric list in buildSeries()
  // (arrays, one entry per row: PLOTLY-STYLE PER-ROW SUBPLOTS), since neither the row
  // count nor titles can live in the JSON-safe serialized skeleton, and since the live
  // count can differ from what Python's SSR skeleton has. One toolbox dataZoom feature
  // spanning every row's xAxisIndex plus one `inside` dataZoom per row (POLISH.md #8):
  // a drag-select inside any row's grid zooms only THAT row's x-axis; double-click reset
  // (wired once per chart in echarts-plotting.js's renderPlot(), generic across every
  // plot type) resets every dataZoom component, all rows included.
  applyOptionOverrides(option) {
    option.grid = this._grids || [];
    option.xAxis = this._xAxis || [];
    option.yAxis = this._yAxis || [];
    const n = this._rowCount || 0;
    option.toolbox = {
      show: true,
      top: "150%",
      feature: {
        dataZoom: { show: true, xAxisIndex: Array.from({ length: n }, (_, i) => i), yAxisIndex: [] },
      },
    };
    option.dataZoom = Array.from({ length: n }, (_, i) => ({
      type: "inside",
      xAxisIndex: [i],
      zoomOnMouseWheel: false,
      moveOnMouseWheel: false,
    }));

    // FIX #3: one combined tooltip per hovered sample, listing that sample's value in
    // EVERY metric row at once (Plotly-style shared hover), not just the hovered row's
    // own value. `_sampleTooltipData` (built in buildSeries() above, alongside each
    // row's scatterData) already holds every row's {title, value, suffix} for a given
    // sample name, so no re-walking of series data is needed here. Falls back to the
    // single hovered-point value on the off chance a sample somehow isn't in the map
    // (e.g. mid-render race), so a tooltip still shows something rather than nothing.
    if (option.tooltip) {
      const self = this;
      option.tooltip.formatter = (params) => {
        if (params.seriesType !== "scatter") return "";
        const rows = self._sampleTooltipData && self._sampleTooltipData[params.name];
        if (!rows || rows.length === 0) {
          let real = Array.isArray(params.value) ? params.value[0] : params.value;
          return `<b>${params.name}</b>: ${window.formatNumber(real)}`;
        }
        const lines = rows.map((r) => `${r.title}: ${window.formatNumber(r.value)}${r.suffix}`);
        return `<b>${params.name}</b><br>${lines.join("<br>")}`;
      };
    }
  }

  // Cross-row sample highlight: hovering a beeswarm point highlights that sample's
  // point on every other metric row. Wired once per chart instance (afterPlotCreated()
  // runs on every render, but the chart instance itself is created once by renderPlot()
  // and reused, so re-registering handlers on every render would stack duplicates).
  afterPlotCreated() {
    if (this._highlightWired || !this.chart) return;
    this._highlightWired = true;
    const chart = this.chart;
    const self = this;

    chart.on("mouseover", { seriesType: "scatter" }, (params) => {
      const targets = self._sampleIndexMap && self._sampleIndexMap[params.name];
      if (!targets || targets.length === 0) return;
      chart.dispatchAction({
        type: "highlight",
        batch: targets.map(([seriesIndex, dataIndex]) => ({ seriesIndex, dataIndex })),
      });
      chart.dispatchAction({ type: "showTip", seriesIndex: params.seriesIndex, dataIndex: params.dataIndex });
    });

    chart.on("mouseout", { seriesType: "scatter" }, (params) => {
      const targets = self._sampleIndexMap && self._sampleIndexMap[params.name];
      if (!targets || targets.length === 0) return;
      chart.dispatchAction({
        type: "downplay",
        batch: targets.map(([seriesIndex, dataIndex]) => ({ seriesIndex, dataIndex })),
      });
    });
  }
}

// Make EchartsViolinPlot globally available (referenced by bare window.EchartsViolinPlot
// in echarts-plotting.js's initPlot(), not a bare identifier, so it's robust across the
// minified bundle regardless of module ordering after bundling).
window.EchartsViolinPlot = EchartsViolinPlot;
