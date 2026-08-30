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
// - `grid`/`xAxis`/`yAxis`/`title` are rebuilt wholesale every render (buildSeries(),
//   applied in applyOptionOverrides()) from the LIVE metric list, since it can differ
//   from the Python skeleton's SSR-time list (e.g. a table column toggled after page
//   load); each row gets a real, SI-formatted x-axis (POLISH.md #6) and a per-row
//   `inside`/toolbox dataZoom so click+drag zoom only affects that row (POLISH.md #8);
//   the row's metric title is drawn by a native ECharts `title` component entry (one per
//   row, appended to `option.title`, see makeRowTitle()), fixed at a pixel position
//   (immune to zoom) and vertically centered (POLISH.md #1) in the left gutter (ROW
//   ALIGNMENT fix: every row's grid now shares the same `left`/`right`/
//   `containLabel: false`, so every row's rect is identical, see gridLeft() below). A
//   `custom` series `renderItem` was tried first (reading `params.coordSys` to stay
//   zoom-independent) but ECharts' SSR export (`static_export.py`) silently drops any
//   manually-drawn `type: "text"` element -- from a `custom` series OR the `graphic`
//   component -- since it has no real Canvas context to measure arbitrary text with;
//   only text from built-in components (title/legend/axis labels) survives SSR.
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
// Measured directly from Plotly's own rendered general-stats violin
// (gsDiv.data[i].fillcolor, live browser inspection): light-mode fillcolor is
// rgba(<color>, 0.5) (templates/default/src/js/plots/violin.js's own
// isDarkMode ? 0.3 : 0.5, and the report always renders in light mode there).
const FILL_ALPHA = 0.5;

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

// Mirrors violin.py::_row_center_pct: vertical center of grid row rowIdx of n, as a
// float percentage of the container height (same basis as rowGeometry()'s top/height,
// duplicated here rather than parsed back out of those formatted strings, same
// cross-language-duplication style as kde()). Used to vertically center that row's
// title (makeRowTitle() below, POLISH.md #1).
function rowCenterPct(rowIdx, n) {
  const total = ROW_PX * n + EXTRA_PX + BOTTOM_PX;
  const topPct = (EXTRA_PX / total) * 100;
  const bottomPct = (BOTTOM_PX / total) * 100;
  const rowSlotPct = (100 - topPct - bottomPct) / n;
  const rowTopPct = topPct + rowIdx * rowSlotPct;
  const rowHeightPct = rowSlotPct * ROW_GRID_FRACTION;
  return rowTopPct + rowHeightPct / 2;
}

// Mirrors violin.py::_grid_left: shared left inset (px) for every row's grid, sized from
// the LONGEST visible title so every row's plot area starts at the same x, computed from
// the title's own character count (not measured) so it never depends on how any one
// row's own label happens to render.
const TITLE_CHAR_PX = 6.0;
const MIN_GRID_LEFT = 60;
const MAX_GRID_LEFT = 280;
// Gap (px) between the right-aligned title text and the row's plot area (ROW ALIGNMENT
// fix, mirrors violin.py::_TITLE_GUTTER_PAD): wide enough to clear the row's OWN x-axis
// first tick label (e.g. "0%"/"0 M"/"627 bp"), which straddles that same left edge (a
// value axis has no boundaryGap) and would otherwise overlap the title text now that
// containLabel/outerBoundsMode no longer reserve room for it (measured worst case ~30px
// wide at fontSize 10, so half its width plus a little slack comfortably fits in 24px).
// See makeRowTitle() below.
const TITLE_GUTTER_PAD = 24;
// Small slack (px) between the title text's own estimated width and where
// TITLE_GUTTER_PAD starts, absorbing TITLE_CHAR_PX's per-glyph estimation error
// (mirrors violin.py::_TITLE_TEXT_SLACK).
const TITLE_TEXT_SLACK = 4;

function gridLeft(titles) {
  if (titles.length === 0) return MIN_GRID_LEFT;
  const maxLen = arrMax(titles.map((t) => t.length));
  const reserved = maxLen * TITLE_CHAR_PX + TITLE_TEXT_SLACK + TITLE_GUTTER_PAD;
  return Math.min(MAX_GRID_LEFT, Math.max(MIN_GRID_LEFT, reserved));
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
// violin.py's _DEFAULT_STROKE/_DEFAULT_FILL. Only used when a metric has no configured
// header.color; a real per-metric color always takes priority below.
const DEFAULT_STROKE = "rgb(153,153,153)";

// Mirrors violin.py::_violin_colors: [fill, stroke], falling back to the same defaults
// (kept as rgb(...), not hex, so toRgba() above applies uniformly). `stroke` is used only
// for the median line (see makeViolinRenderItem): Plotly's own violin measures
// line.width: 0 for its outline/box/meanline (verified via its rendered _fullData), so
// neither the polygon nor the box gets a border here.
function violinColors(header) {
  const rgb = header.color ? normalizeColorToRGB(header.color) : null;
  if (!rgb) return [toRgba(DEFAULT_STROKE, FILL_ALPHA), DEFAULT_STROKE];
  const stroke = `rgb(${rgb})`;
  return [toRgba(stroke, FILL_ALPHA), stroke];
}

// Draws the KDE polygon plus an inner Q1-Q3 box and median line (POLISH.md #10b),
// matching Plotly's box_visible/meanline violin config. q1/median/q3 are real values
// (same x-space as poly). Mirrors violin.py::_violin_render_series: neither the polygon
// nor the box has a stroke (matches Plotly's measured line.width: 0); the box reuses the
// SAME fill as the body (also matching Plotly, whose box fillcolor equals the trace's
// own), reading as a slightly denser patch purely from the two semi-transparent fills
// overlapping. Only the median line keeps a solid stroke (Plotly's own meanline.width
// measures 0, i.e. invisible; a fully invisible median would be a regression, not a
// match, so it's kept here, just thin).
function makeViolinRenderItem(poly, q1, median, q3, fill, stroke) {
  return function (params, api) {
    const pts = poly.map((p) => api.coord(p));
    const bTL = api.coord([q1, -BOX_HALF_HEIGHT]);
    const bBR = api.coord([q3, BOX_HALF_HEIGHT]);
    const mTop = api.coord([median, -MEDIAN_HALF_HEIGHT]);
    const mBot = api.coord([median, MEDIAN_HALF_HEIGHT]);
    return {
      type: "group",
      children: [
        { type: "polygon", shape: { points: pts }, style: { fill } },
        {
          type: "rect",
          shape: {
            x: Math.min(bTL[0], bBR[0]),
            y: Math.min(bTL[1], bBR[1]),
            width: Math.abs(bBR[0] - bTL[0]),
            height: Math.abs(bBR[1] - bTL[1]),
          },
          style: { fill },
        },
        {
          type: "line",
          shape: { x1: mTop[0], y1: mTop[1], x2: mBot[0], y2: mBot[1] },
          style: { stroke, lineWidth: 2 },
        },
      ],
    };
  };
}

// One entry in ECharts' native `title` component array (option.title) drawing the row's
// metric TITLE, fixed at gridLeft - TITLE_GUTTER_PAD px from the container's left edge
// and vertically centered on the row (rowCenterPct), right-aligned so its text ends
// exactly there. padding:0/fontWeight:"normal" strip the title component's own defaults
// (a title normally reserves internal padding and renders bold).
//
// A NATIVE ECharts component, not a `custom` series renderItem or `graphic` element:
// tried both of those first, since either can read a fixed pixel position immune to zoom
// (see this file's header comment), but ECharts' SSR mode (static_export.py) silently
// drops manually-drawn type:"text" elements from EITHER of those, because it has no real
// Canvas 2D context to measure arbitrary text with. `title` array entries are a built-in
// component (same family as the chart's own main title, axis labels, and legend) that
// ships its own SSR-safe measurement fallback, so this is the one placement verified to
// render in both the browser and the MiniRacer/resvg SSR path. Mirrors
// violin.py::_title_option.
function makeRowTitle(rowIdx, n, title, gridLeftPx, textColor) {
  return {
    text: title,
    left: gridLeftPx - TITLE_GUTTER_PAD,
    top: `${rowCenterPct(rowIdx, n)}%`,
    textAlign: "right",
    textVerticalAlign: "middle",
    padding: 0,
    textStyle: { fontSize: 12, fontWeight: "normal", color: textColor },
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
  // n_metric beeswarm `scatter` series -- same ordering as the SSR path
  // (multiqc/plots/echarts/violin.py::series()), which is what the sample ->
  // [seriesIndex, dataIndex] cross-highlight map below relies on. There is no min/max
  // text-annotation series any more (removed: each row already has its own real-valued
  // x-axis, so a text label repeating the same range read as if the beeswarm dots
  // themselves carried on-canvas labels). Also (re)builds
  // `this._grids`/`this._xAxis`/`this._yAxis`/`this._titles`/
  // `this._toolbox`/`this._dataZoom` from the LIVE metric list (PLOTLY-STYLE PER-ROW
  // SUBPLOTS), applied onto the option in applyOptionOverrides() below. Row titles
  // (`this._titles`) are native `title` component entries, not a series (see
  // makeRowTitle()).
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
      this._titles = [];
      this._rowCount = 0;
      return [];
    }

    const colors = window.getEchartsThemeColors();
    const rowLeft = gridLeft(rows.map(({ header }) => metricTitle(header)));
    const violinSeries = [];
    const titles = [];
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

      const [fill, stroke] = violinColors(header);
      violinSeries.push({
        type: "custom",
        coordinateSystem: "cartesian2d",
        xAxisIndex: rowIdx,
        yAxisIndex: rowIdx,
        data: [0],
        renderItem: makeViolinRenderItem(poly, q1, median, q3, fill, stroke),
        silent: true,
        z: 1,
      });

      const title = metricTitle(header);
      titles.push(makeRowTitle(rowIdx, n, title, rowLeft, colors.tickcolor));

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
        // Above the violin body's z: 1 (KDE polygon/box/median), so points always draw
        // on top of the (now un-bordered, semi-transparent) fill, never obscured by it.
        z: 3,
      });

      // Grid/axis geometry for this row (PLOTLY-STYLE PER-ROW SUBPLOTS): a real,
      // SI-formatted x-axis (POLISH.md #6), and a FULLY HIDDEN value y-axis (just a
      // [-0.5, 0.5] coordinate range for the renderItems above to draw the violin's
      // thickness against; the row title is drawn separately by makeRowTitle(), not by
      // this axis, see the ROW ALIGNMENT fix comment on gridLeft() above).
      // `left`/`right` are IDENTICAL for every row and `containLabel: false` means
      // nothing can grow one row's rect differently from another's, so every row's
      // coordinateSystem.getRect() ends up equal in both x and width. Axis colors are
      // themed to match every other plot type (mirrors the buildCurrentOption theme step
      // in echarts-plotting.js, which never runs on these axes since they're built fresh
      // here, after that step already ran).
      const [top, gridHeight] = rowGeometry(rowIdx, n);
      grids.push({
        top,
        height: gridHeight,
        left: rowLeft,
        right: 16,
        containLabel: false,
        // ECharts 6 also auto-nudges a grid's rect inward, independently of
        // containLabel, to keep axis labels from overflowing the chart's OUTER bounds
        // (grid.outerBoundsMode default "auto", see GridModel.js's OUTER_BOUNDS_DEFAULT):
        // since a row's x-axis tick label LENGTH varies by metric (e.g. "20%" vs
        // "18.39 M"), this on its own reproduces the exact same per-row rect inequality
        // bug containLabel caused, just via a different ECharts 6 feature. "none"
        // disables it, verified empirically to make every row's rect byte-for-byte
        // identical (mirrors violin.py::_row_axes).
        outerBoundsMode: "none",
      });
      xAxis.push({
        type: "value",
        gridIndex: rowIdx,
        min: lo,
        max: hi,
        axisLabel: { fontSize: 10, color: colors.tickcolor, formatter: (v) => window.formatAxisNumber(v, suffix) },
        // No axis line/ticks (measured Plotly's own violin x-axis via its _fullLayout:
        // showline: false, ticks: ""; only the tick label text is shown). A full-width
        // axis line plus tick marks under every row read far more present/"brighter"
        // than Plotly's bare label text, and was the main source of the axis-brightness
        // mismatch.
        axisLine: { show: false },
        axisTick: { show: false },
        splitLine: { show: false },
      });
      yAxis.push({
        type: "value",
        gridIndex: rowIdx,
        min: -0.5,
        max: 0.5,
        interval: 0.5,
        axisLabel: { show: false },
        axisTick: { show: false },
        axisLine: { show: false },
        splitLine: { show: false },
      });
    });

    this._grids = grids;
    this._xAxis = xAxis;
    this._yAxis = yAxis;
    this._titles = titles;
    this._rowCount = n;

    // sample -> [[seriesIndex, dataIndex], ...] across every scatter series, for the
    // cross-row highlight wired in afterPlotCreated().
    const sampleIndexMap = {};
    const scatterOffset = violinSeries.length;
    scatterSeries.forEach((series, si) => {
      series.data.forEach((item, di) => {
        if (!sampleIndexMap[item.name]) sampleIndexMap[item.name] = [];
        sampleIndexMap[item.name].push([scatterOffset + si, di]);
      });
    });
    this._sampleIndexMap = sampleIndexMap;
    this._sampleTooltipData = sampleTooltipData;

    return violinSeries.concat(scatterSeries);
  }

  // grid/xAxis/yAxis/title are rebuilt wholesale from the live metric list in
  // buildSeries() (arrays, one entry per row: PLOTLY-STYLE PER-ROW SUBPLOTS), since
  // neither the row count nor titles can live in the JSON-safe serialized skeleton, and
  // since the live count can differ from what Python's SSR skeleton has. Row titles
  // (this._titles) are appended after the chart's own main title (already themed by
  // buildCurrentOption's generic step, which runs before this method, see
  // echarts-plotting.js), not replacing it: ECharts' `title` option accepts either one
  // component or a list of them. One toolbox dataZoom feature spanning every row's
  // xAxisIndex plus one `inside` dataZoom per row (POLISH.md #8): a drag-select inside
  // any row's grid zooms only THAT row's x-axis; double-click reset (wired once per chart
  // in echarts-plotting.js's renderPlot(), generic across every plot type) resets every
  // dataZoom component, all rows included.
  applyOptionOverrides(option) {
    option.grid = this._grids || [];
    option.xAxis = this._xAxis || [];
    option.yAxis = this._yAxis || [];
    option.title = (option.title ? [option.title] : []).concat(this._titles || []);
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
      // filterMode: "none" (default is "filter"): ECharts' default dataZoom behavior
      // REMOVES a series' data points once their x-value falls outside the zoomed range,
      // which for the violin custom series (a dummy data: [0], not a real x-coordinate;
      // its renderItem draws from closure state, not from that data value) makes ECharts
      // drop the whole row's polygon the moment a zoom narrows the axis so 0 falls
      // outside it. "none" zooms the axis view without filtering any series data
      // (mirrors violin.py::layout_option).
      filterMode: "none",
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
