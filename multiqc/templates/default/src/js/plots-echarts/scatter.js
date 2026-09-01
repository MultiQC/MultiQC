// Plotly marker symbol names (as used by point.marker_symbol) mapped to their closest
// ECharts built-in symbol. ECharts' built-in set is much smaller than Plotly's (no
// per-symbol open/dot variants, no exotic shapes), so this only covers the common names
// with a clean visual equivalent; anything else falls back to "circle", same as the
// Python side (multiqc/plots/echarts/scatter.py's `_echarts_symbol()`).
const PLOTLY_TO_ECHARTS_SYMBOL = {
  circle: "circle",
  square: "rect",
  diamond: "diamond",
  triangle: "triangle",
  "triangle-up": "triangle",
};

function echartsSymbol(plotlySymbol) {
  if (!plotlySymbol) return "circle";
  return PLOTLY_TO_ECHARTS_SYMBOL[plotlySymbol] ?? "circle";
}

// ECharts scatter plot: extends the default template's ScatterPlot (imported just before
// this file in main-js.js) purely for its prepData()/exportData()/formatDatasetForAiPrompt()/
// activeDatasetSize(), which are engine-neutral (they only read/filter raw dataset fields
// through applyToolboxSettings()). window.ScatterPlot itself extends window.Plot, which by
// this point in the import order is OUR echarts Plot class (from echarts-plotting.js), not
// Plotly's.
//
// Series strategy: when no point carries a `group` field, ONE `{type: "scatter"}` series
// per dataset with per-item styling, not one series per point (see
// multiqc/plots/echarts/scatter.py docstring for why: it keeps a large point cloud a
// single series, which is what makes the canvas-vs-svg renderer threshold meaningful).
// When points DO carry a `group` field, one series PER GROUP instead (see
// buildGroupedSeries), so ECharts' native per-series legend gives Plotly's
// legendgroup/inLegend behaviour: one legend entry per group, click-to-toggle the whole
// group. Real grouped datasets have few groups, so this doesn't defeat the threshold.
class EchartsScatterPlot extends window.ScatterPlot {
  buildSeries() {
    let dataset = this.datasets[this.activeDatasetIdx];
    let marker = dataset["trace_params"]?.["marker"] || {};
    let defaultColor = marker.color ?? "rgba(124, 181, 236, .5)";
    let defaultSize = marker.size ?? 10;

    let [samples, points] = this.prepData();
    // prepData() (inherited from default ScatterPlot) leaves a hole (undefined) for
    // every point the toolbox hides instead of removing it. Compute the "N samples
    // hidden" banner from the hole count before dropping them, same pattern as bar.js.
    this._updateHiddenSamplesWarning(samples.length, points.filter(Boolean).length);
    points = points.filter(Boolean);

    if (points.length === 0 || samples.length === 0) return [];

    // Reorder so highlighted points draw on top, same as default plots/scatter.js.
    let highlighted = points.filter((p) => p.highlight);
    let nonHighlighted = points.filter((p) => !p.highlight);
    points = nonHighlighted.concat(highlighted);

    // Sort points by group order if groups are provided in config, same as default
    // plots/scatter.js buildTraces(). This also determines the group series order below
    // (first-seen order of a now-group-sorted point list matches pconfig.groups).
    if (this.pconfig.groups) {
      points.sort((a, b) => {
        let aIndex = this.pconfig.groups.indexOf(a.group);
        let bIndex = this.pconfig.groups.indexOf(b.group);
        if (aIndex === -1) aIndex = this.pconfig.groups.length;
        if (bIndex === -1) bIndex = this.pconfig.groups.length;
        return aIndex - bIndex;
      });
    }

    // Highlight rule keyed off the GLOBAL flag (not just this plot's matches), matching
    // the dim/grey convention landed for other plot types: once any highlight is active
    // anywhere, every non-matching point greys out.
    let anyHighlighted = window.mqc_highlight_f_texts.length > 0;
    let buildItem = (point) => {
      let color = anyHighlighted ? (point.highlight ?? "#cccccc") : (point.color ?? defaultColor);
      let itemStyle = { color: color };
      if (point.opacity !== undefined) itemStyle.opacity = point.opacity;
      if (point.marker_line_width !== undefined) itemStyle.borderWidth = point.marker_line_width;

      return {
        value: [point.x, point.y],
        name: point.name,
        symbolSize: point.marker_size ?? defaultSize,
        symbol: echartsSymbol(point.marker_symbol),
        itemStyle: itemStyle,
      };
    };

    let hasGroups = points.some((point) => point.group);
    let series = hasGroups
      ? this.buildGroupedSeries(points, buildItem)
      : [
          {
            type: "scatter",
            name: dataset["label"] ?? "",
            data: points.map(buildItem),
          },
        ];

    series.push(...this.bandsLinesSeries());
    return series;
  }

  // One series per group, `name` = group name, so ECharts' legend shows one entry per
  // group and toggling it hides/shows the whole group. Points are bucketed by
  // point.group; a point with no group falls back to its own point.name as the bucket
  // key, the same fallback the default template's Plotly buildTraces() inLegend key
  // uses, so a point with no group still gets its own legend entry rather than being
  // silently folded into another one. `points` is expected to already be sorted by
  // pconfig.groups (buildSeries() does this above), so the first-seen bucket order here
  // already matches.
  buildGroupedSeries(points, buildItem) {
    let order = [];
    let buckets = new Map();
    points.forEach((point) => {
      let key = point.group ?? point.name;
      if (!buckets.has(key)) {
        buckets.set(key, []);
        order.push(key);
      }
      buckets.get(key).push(point);
    });

    return order.map((key) => ({
      type: "scatter",
      name: key,
      data: buckets.get(key).map(buildItem),
    }));
  }

  // See EchartsBarPlot._updateHiddenSamplesWarning (plots-echarts/bar.js) for the full
  // rationale; copied here (kept minimal) since ScatterPlot and BarPlot don't share a
  // common ECharts base beyond window.Plot.
  _updateHiddenSamplesWarning(total, visible) {
    let groupDiv = $("#" + this.anchor).closest(".mqc_hcplot_plotgroup");
    groupDiv.parent().find(".samples-hidden-warning").remove();

    if (visible === 0 && total > 0) {
      groupDiv.hide();
      return;
    }
    groupDiv.show();

    let nHidden = total - visible;
    if (nHidden > 0) {
      const alert = `
      <div class="samples-hidden-warning alert alert-warning">
        ⚠ <strong>Warning:</strong> ${nHidden} samples hidden.
        <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose('#mqc_hidesamples', true); return false;">See toolbox.</a>
      </div>`;
      groupDiv.before(alert);
    }
  }

  // Tooltip formatter can't live in the JSON-safe serialized skeleton, so it's attached
  // here to the fully-assembled option (see Plot.applyOptionOverrides in
  // echarts-plotting.js).
  applyOptionOverrides(option) {
    if (!option.tooltip) return;
    let xsuffix = this.layout?.xaxis?.ticksuffix ?? "";
    let ysuffix = this.layout?.yaxis?.ticksuffix ?? "";
    // Matches the Plotly reference's native hovertemplate, which applies
    // `layout.xaxis.hoverformat`/`layout.yaxis.hoverformat` to `%{x}`/`%{y}` respectively.
    // NOTE: read off the active dataset's raw layout dict, not `this.layout` (see
    // bar.js's matching comment): the neutral `AxisIR` behind `this.layout` drops
    // `hoverformat`, a Plotly-only key.
    let dsLayout = this.datasets[this.activeDatasetIdx]?.layout;
    let xFormat = dsLayout?.xaxis?.hoverformat;
    let yFormat = dsLayout?.yaxis?.hoverformat;
    option.tooltip.formatter = (params) => {
      let [x, y] = params.value;
      return `<b>${params.name}</b><br/>X: ${window.formatWithHoverformat(x, xFormat)}${xsuffix}<br/>Y: ${window.formatWithHoverformat(y, yFormat)}${ysuffix}`;
    };
  }
}

// Make EchartsScatterPlot globally available (referenced by bare window.EchartsScatterPlot
// in echarts-plotting.js's initPlot(), not a bare identifier, so it's robust across the
// minified bundle regardless of module ordering after bundling).
window.EchartsScatterPlot = EchartsScatterPlot;
