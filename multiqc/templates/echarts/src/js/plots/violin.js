// ECharts violin plot: STANDALONE class extending window.Plot directly (Task 2.2).
// Unlike bar/line/scatter/box, this does NOT extend the default template's ViolinPlot:
// that class's buildTraces() builds one Plotly subplot (its own x/y axis pair) per
// metric row, which has no ECharts equivalent (a single value-axis-trick grid is used
// instead, see multiqc/plots/echarts/violin.py's module docstring). Only the
// engine-neutral bits are ported here verbatim: `prepData()` field access
// (templates/default/src/js/plots/violin.js:18-43) and `exportData()` (:401-431).
//
// Rendering strategy (mirrors multiqc/plots/echarts/violin.py):
// - Each row is a `custom` series whose renderItem draws a closed KDE polygon.
// - `plot.echarts.datasets[i].violins` ({metric: {poly, range}}) is precomputed by
//   Python from ALL samples. When no sample is hidden by the toolbox we reuse that
//   polygon directly (translating its y-offset if the row position shifted, e.g. a
//   metric got hidden via the table column config); when a sample IS hidden we
//   recompute the KDE in JS from the visible values only (kde() below, a byte-for-byte
//   port of violin.py::kde -- see the golden fixture comment above it).
// - A beeswarm `scatter` series per row uses ECharts 6's native `jitter`.
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
const SCATTER_COLOR = "rgba(30,50,80,0.85)";

// Reduce-based min/max: `Math.min(...arr)`/`Math.max(...arr)` blow the call stack on
// large sample counts (general stats tables can have thousands of rows).
function arrMin(arr) {
  return arr.reduce((a, b) => (b < a ? b : a), arr[0]);
}
function arrMax(arr) {
  return arr.reduce((a, b) => (b > a ? b : a), arr[0]);
}

// Mirrors violin.py::_metric_range (RANGE_PAD-padded [lo, hi] domain).
function metricRange(values) {
  let lo = arrMin(values);
  let hi = arrMax(values);
  const span = hi - lo;
  lo -= span * RANGE_PAD;
  hi += span * RANGE_PAD;
  if (hi <= lo) hi = lo + 1.0;
  return [lo, hi];
}

// Mirrors violin.py::_violin_polygon: closed KDE polygon in normalized (0..1) x-space
// and integer-offset y-space.
function violinPolygon(values, lo, hi, rowIdx) {
  const span = hi - lo;
  const xs = Array.from({ length: N_BINS }, (_, i) => lo + (span * i) / (N_BINS - 1));
  const ys = kde(values, xs);
  const ymax = arrMax(ys) || 1.0;
  const top = xs.map((x, i) => [(x - lo) / span, rowIdx + (ys[i] / ymax) * ROW_HEIGHT]);
  const bottom = xs.map((x, i) => [(x - lo) / span, rowIdx - (ys[i] / ymax) * ROW_HEIGHT]).reverse();
  return top.concat(bottom);
}

// Mirrors violin.py::_metric_title.
function metricTitle(header) {
  return header.namespace ? `${header.namespace}: ${header.title}` : header.title;
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

// Mirrors violin.py::_violin_colors: [fill, stroke], falling back to the same defaults.
function violinColors(header) {
  const rgb = header.color ? normalizeColorToRGB(header.color) : null;
  if (!rgb) return ["rgba(91,143,249,0.5)", "#1f4d9e"];
  return [`rgba(${rgb},0.5)`, `rgb(${rgb})`];
}

// Min/max row label: delegates to the shared window.formatNumber (echarts-plotting.js)
// so every echarts plot type rounds numbers the same way.
function formatG(num) {
  return String(window.formatNumber(num));
}

function makeKdeRenderItem(poly, fill, stroke) {
  return function (params, api) {
    const pts = poly.map((p) => api.coord(p));
    return { type: "polygon", shape: { points: pts }, style: { fill, stroke, lineWidth: 1 } };
  };
}

function makeAnnotationRenderItem(rowIdx, loObs, hiObs) {
  const loLabel = formatG(loObs);
  const hiLabel = formatG(hiObs);
  return function (params, api) {
    const L = api.coord([0, rowIdx]);
    const R = api.coord([1, rowIdx]);
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
    this.violinHeight = dump["violin_height"];
    this.extraHeight = dump["extra_height"];
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

  // Builds: n_metric KDE-polygon `custom` series, then n_metric min/max-annotation
  // `custom` series, then n_metric beeswarm `scatter` series -- same ordering as the
  // SSR path (multiqc/plots/echarts/violin.py::series()), which is what the
  // sample -> [seriesIndex, dataIndex] cross-highlight map below relies on.
  buildSeries() {
    let [
      metrics,
      headerByMetric,
      allSamples,
      sampleSettings,
      violinValuesBySampleByMetric,
      scatterValuesBySampleByMetric,
    ] = this.prepData();
    if (metrics.length === 0) return [];

    const someHidden = sampleSettings.some((s) => s.hidden);
    const highlightingEnabled = sampleSettings.some((s) => s.highlight !== null);

    const echartsDs = this.echarts.datasets[this.activeDatasetIdx];
    const violinsFromPython = echartsDs.violins || {};
    // Row index each metric occupied when Python precomputed the polygons (from ALL
    // samples, keyed by insertion order): needed to re-align a reused polygon's baked-in
    // y-offset when the CURRENT row position differs (e.g. a metric got hidden via the
    // table column config, which Python's serialization doesn't know about).
    const pythonRowIdxByMetric = {};
    Object.keys(violinsFromPython).forEach((m, i) => {
      pythonRowIdxByMetric[m] = i;
    });

    const kdeSeries = [];
    const annotationSeries = [];
    const scatterSeries = [];
    const rowTitles = {};
    let rowIdx = 0;

    metrics.forEach((metric) => {
      const header = headerByMetric[metric];
      const valuesBySample = violinValuesBySampleByMetric[metric] || {};
      const numericValues = Object.values(valuesBySample).filter((v) => typeof v === "number" && Number.isFinite(v));
      if (numericValues.length === 0) return; // no data left to draw this row (e.g. all its samples hidden)

      let poly, lo, hi;
      const pre = violinsFromPython[metric];
      if (!someHidden && pre) {
        // Fast path: reuse the precomputed polygon, translating its y-offset if this
        // metric's row position shifted relative to Python's serialization order.
        [lo, hi] = pre.range;
        const pythonRowIdx = pythonRowIdxByMetric[metric] ?? rowIdx;
        const dy = rowIdx - pythonRowIdx;
        poly = dy === 0 ? pre.poly : pre.poly.map(([x, y]) => [x, y + dy]);
      } else {
        // Slow path: samples are hidden, so the polygon must be recomputed from the
        // currently-visible values only (see the perf note in BUILD_PLAN.md Phase 2 risks).
        [lo, hi] = metricRange(numericValues);
        poly = violinPolygon(numericValues, lo, hi, rowIdx);
      }

      const [fill, stroke] = violinColors(header);
      kdeSeries.push({
        type: "custom",
        coordinateSystem: "cartesian2d",
        data: [0],
        renderItem: makeKdeRenderItem(poly, fill, stroke),
        silent: true,
        z: 1,
      });

      annotationSeries.push({
        type: "custom",
        coordinateSystem: "cartesian2d",
        data: [0],
        renderItem: makeAnnotationRenderItem(rowIdx, arrMin(numericValues), arrMax(numericValues)),
        silent: true,
        z: 3,
      });

      rowTitles[rowIdx] = metricTitle(header);

      const scatterData = [];
      if (header.show_points) {
        const svBySample = scatterValuesBySampleByMetric[metric] || {};
        const span = hi - lo;
        Object.entries(svBySample).forEach(([sample, value]) => {
          if (typeof value !== "number" || !Number.isFinite(value)) return;
          const normx = (value - lo) / span;
          const state = sampleSettings[allSamples.indexOf(sample)];
          let color = SCATTER_COLOR;
          let size = 6;
          if (highlightingEnabled) {
            color = state?.highlight ?? "#ddd";
            size = state?.highlight != null ? 8 : 6;
          }
          scatterData.push({
            // Real (denormalized) value stored as the 3rd element for the tooltip.
            value: [normx, rowIdx, value],
            name: state?.name ?? sample,
            itemStyle: { color, borderColor: "#fff", borderWidth: 0.5 },
            symbolSize: size,
          });
        });
      }
      scatterSeries.push({
        type: "scatter",
        name: String(metric),
        data: scatterData,
        symbolSize: 6,
        jitter: 22,
        jitterOverlap: false,
        z: 2,
      });

      rowIdx++;
    });

    // sample -> [[seriesIndex, dataIndex], ...] across every scatter series, for the
    // cross-row highlight wired in afterPlotCreated().
    const sampleIndexMap = {};
    const scatterOffset = kdeSeries.length + annotationSeries.length;
    scatterSeries.forEach((series, si) => {
      series.data.forEach((item, di) => {
        if (!sampleIndexMap[item.name]) sampleIndexMap[item.name] = [];
        sampleIndexMap[item.name].push([scatterOffset + si, di]);
      });
    });
    this._sampleIndexMap = sampleIndexMap;
    this._rowTitles = rowTitles;
    this._rowCount = rowIdx;

    // Wrapper/container height depends on the current row count (metrics can be
    // hidden/shown live via the table column config), same as default violin.js's
    // buildTraces() setting `$(wrapper).css("height", ...)`.
    const height = (this.violinHeight || 70) * Math.max(rowIdx, 1) + (this.extraHeight || 63);
    const el = document.getElementById(this.anchor);
    if (el) el.style.height = height + "px";
    $("#" + this.anchor + "-wrapper").css("height", height + "px");

    return kdeSeries.concat(annotationSeries, scatterSeries);
  }

  // Row labels (yAxis is a value axis, per the VALUE-AXIS TRICK in violin.py) and the
  // tooltip formatter: neither can live in the JSON-safe serialized skeleton.
  applyOptionOverrides(option) {
    const rowTitles = this._rowTitles || {};
    if (option.yAxis) {
      option.yAxis.axisLabel = { ...option.yAxis.axisLabel, formatter: (v) => rowTitles[Math.round(v)] || "" };
      option.yAxis.max = Math.max((this._rowCount || 0) - 0.5, 0.5);
    }
    if (option.tooltip) {
      option.tooltip.formatter = (params) => {
        if (params.seriesType !== "scatter") return "";
        let real = Array.isArray(params.value) ? params.value[2] : params.value;
        return `<b>${params.name}</b>: ${window.formatNumber(real)}`;
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
