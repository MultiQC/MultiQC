// ECharts seqcontent plot: extends the default template's SeqContentPlot (imported just
// before this file in main-js.js) for its prepData()/exportData()/formatDatasetForAiPrompt()/
// activeDatasetSize(), which are engine-neutral (they only read/filter raw dataset fields
// through applyToolboxSettings()). window.SeqContentPlot itself extends window.Plot, which by
// this point in the import order is OUR echarts Plot base class (echarts-plotting.js), not
// Plotly's. Unlike heatmap.js (standalone), the default seqcontent.js has no top-level
// `$(function () {...})` handler calling Plotly directly, so it's safe to import and extend,
// same strategy as bar/line/scatter/box.

// GOLDEN CONTRACT: kept in lockstep with `bin_rgb()` in multiqc/plots/seqcontent.py and
// `seqContentBinRgb()` in templates/default/src/js/plots/seqcontent.js:
// R = %T, G = %A, B = %C (%G implied by the complement of the other three).
function seqContentBinRgb(bin) {
  return [Math.round((bin.t / 100) * 255), Math.round((bin.a / 100) * 255), Math.round((bin.c / 100) * 255)];
}

// GOLDEN CONTRACT: this is the LIVE JS twin of `RENDER_ITEM_BODY` in
// multiqc/plots/echarts/seqcontent.py (the `__FN__` sentinel body executed by
// static_export.py's SSR walker for the flat/kaleido-free export path). Both must map a
// `[start, end, row, r, g, b, opacity]` data item to the same rect geometry, fill color,
// and opacity; a bin covers columns `start..end` inclusive (1-based), so its rect spans
// the data-space interval `[start, end + 1)` on the value xAxis, and `row` (a
// category-axis index) is used as-is. The rect is padded 1px past its true data-space
// edges on both axes (width/height +1, height still centered on the row) so neighbouring
// bins overlap by ~1px instead of leaving a hairline seam, the way the old canvas
// renderer's solid per-pixel fill never showed gaps; no stroke is ever set, so there is
// no border to remove either. `opacity` is the toolbox-highlight dim factor computed by
// buildSeries() below.
//
// CLIPPING: with `filterMode: "none"` (seqcontent.py's `layout_option`), dataZoom only
// moves the axis view range and never removes out-of-range items from `data`, so a bin
// straddling the zoomed edge still reaches renderItem with its untouched geometry. The
// SVG renderer does not auto-clip custom series to the grid the way it clips built-in
// series (bar/line/etc), so an unclipped rect paints straight over the y-axis sample
// labels and the plot title. Rather than depend on a global `echarts.graphic.clipRectByRect`
// (the Python twin's SSR path cannot assume that's in scope: static_export.py executes its
// copy of this body via MiniRacer against a fresh V8 context), the rect is intersected
// against `params.coordSys` (the cartesian grid's pixel rect: `{x, y, width, height}`) by
// hand: shrink to the overlap, or draw nothing (`return;`, i.e. `undefined`, which ECharts
// treats as "no element for this data item") when the rect falls fully outside the grid.
function seqContentRenderItem(params, api) {
  var start = api.value(0);
  var end = api.value(1);
  var row = api.value(2);
  var r = api.value(3);
  var g = api.value(4);
  var b = api.value(5);
  var opacity = api.value(6);
  var p0 = api.coord([start, row]);
  var p1 = api.coord([end + 1, row]);
  var height = api.size([0, 1])[1] + 1;
  var width = p1[0] - p0[0] + 1;
  var rx0 = p0[0];
  var ry0 = p0[1] - height / 2;
  var rx1 = rx0 + width;
  var ry1 = ry0 + height;
  var grid = params.coordSys;
  var gx0 = grid.x;
  var gy0 = grid.y;
  var gx1 = gx0 + grid.width;
  var gy1 = gy0 + grid.height;
  var cx0 = Math.max(rx0, gx0);
  var cy0 = Math.max(ry0, gy0);
  var cx1 = Math.min(rx1, gx1);
  var cy1 = Math.min(ry1, gy1);
  if (cx1 <= cx0 || cy1 <= cy0) return;
  return {
    type: "rect",
    shape: { x: cx0, y: cy0, width: cx1 - cx0, height: cy1 - cy0 },
    style: { fill: "rgb(" + r + "," + g + "," + b + ")", opacity: opacity },
  };
}

class EchartsSeqContentPlot extends window.SeqContentPlot {
  // ECharts analog of default seqcontent.js's buildTraces(): one bin -> one rect, with its
  // true variable width (no expansion into a uniform per-bp pixel grid, unlike the Plotly
  // side; see multiqc/plots/echarts/seqcontent.py's module docstring).
  buildSeries() {
    let [sampleSettings, rows] = this.prepData();
    if (sampleSettings.length === 0) return [];

    // Stashed for buildCurrentOption() (yAxis, generic single-axis path) and
    // applyOptionOverrides() below (tooltip needs the sample names and bins too).
    this._axisData = { axis: "yAxis", data: sampleSettings.map((s) => s.name) };
    this._rows = rows;

    // Highlight (POLISH item 6): when any sample is highlighted via the toolbox, dim
    // every non-highlighted row's rects instead of leaving the plot unchanged, same
    // dim-non-highlighted rule as bar.js's alpha handling. Re-runs on every render, so
    // it stays current with the toolbox's `mqc_highlights` event (plotting-shared.js
    // already calls window.renderPlot() for every plot on that event).
    let highlighted = sampleSettings.filter((s) => s.highlight);
    let data = [];
    for (let rowIdx = 0; rowIdx < rows.length; rowIdx++) {
      let opacity = highlighted.length > 0 && !sampleSettings[rowIdx].highlight ? 0.25 : 1;
      for (let bin of rows[rowIdx]) {
        let [r, g, b] = seqContentBinRgb(bin);
        data.push([bin.start, bin.end, rowIdx, r, g, b, opacity]);
      }
    }

    return [
      {
        type: "custom",
        coordinateSystem: "cartesian2d",
        renderItem: seqContentRenderItem,
        data: data,
        // REGRESSION FIX (Item 1): see the matching comment on the `encode` field of
        // the SSR twin's series() in multiqc/plots/echarts/seqcontent.py -- without it,
        // ECharts infers dim0 -> x, dim1 -> y for a cartesian custom series, so the
        // `zoom_option` yAxis dataZoom filtered items by `end` (item index 1) against
        // the sample-count range instead of by `row` (item index 2), silently dropping
        // every bin whose `end` exceeded the sample count. This encode fixes that.
        encode: { x: [0, 1], y: 2 },
        // Item 4: pointer cursor over this heatmap's rects (its primary interaction is
        // click-to-drill-down, see afterPlotCreated() below), not ECharts' default
        // zoom/crosshair cursor. `cursor` is a series-level option, so it only affects
        // this series, never other echarts plot types.
        cursor: "pointer",
      },
    ];
  }

  // Tooltip formatter can't live in the JSON-safe serialized skeleton, so it's attached
  // here to the fully-assembled option (see Plot.applyOptionOverrides in
  // echarts-plotting.js). Looks the bin back up by its data-item start/end rather than
  // reading params.data[3..5] (the precomputed RGB) so the readout always shows the
  // original a/c/g/t percentages, not a value derived from the color.
  applyOptionOverrides(option) {
    if (!option.tooltip) return;
    let axisData = this._axisData ? this._axisData.data : [];
    let rows = this._rows ?? [];
    // Colored swatch before each base's readout (POLISH item 2), same T/C/A/G ->
    // color mapping the Plotly template's tooltip uses.
    let swatch = (color) =>
      `<span style="display:inline-block;width:10px;height:10px;background:${color};margin-right:4px;"></span>`;
    option.tooltip.formatter = (params) => {
      let [start, end, rowIdx] = params.value;
      let bin = (rows[rowIdx] ?? []).find((b) => b.start === start && b.end === end);
      if (!bin) return "";
      let sampleName = axisData[rowIdx] ?? "";
      let posLabel = start === end ? `${start} bp` : `${start}-${end} bp`;
      return (
        `<b>${sampleName}</b><br/>${posLabel}<br/>` +
        `${swatch("#dc0000")}%T: ${bin.t}%<br/>${swatch("#0000dc")}%C: ${bin.c}%<br/>` +
        `${swatch("#00dc00")}%A: ${bin.a}%<br/>${swatch("#404040")}%G: ${bin.g}%`
      );
    };
  }

  // T3.1 click-to-drilldown: bind ECharts' click event to the SAME openDrilldown()/
  // bindDrilldownControls() shared logic defined on window.SeqContentPlot (the default
  // template's class, our superclass here) -- no swap/show/render logic duplicated.
  // Custom-series click data items are `[start, end, row, r, g, b, opacity]`
  // (buildSeries() above), so params.data[2] is the row index into the SAME toolbox-filtered rows
  // prepData() computed for this render.
  afterPlotCreated() {
    this.bindDrilldownControls();
    if (!this._cursorStyleInjected) {
      this._cursorStyleInjected = true;
      // Item 4: the series-level `cursor: "pointer"` set in buildSeries() above only
      // covers the rendered rects themselves; ECharts' box-zoom "dataZoomSelect" cursor
      // mode (activated globally by echarts-plotting.js's renderPlot() whenever
      // `toolbox.feature.dataZoom` is present, which it is here, see zoom_option) sets
      // an inline `cursor: crosshair` on the chart's own wrapping div on every render,
      // which otherwise wins over empty space between rects (and, since it's inline,
      // over any plain CSS rule without `!important`). Scoped by this plot's own
      // anchor id (never touches other echarts plots), `!important` beats that inline
      // style the same way the default template's Item 4 fix beats Plotly's
      // `.nsewdrag` inline cursor.
      const style = document.createElement("style");
      style.textContent = `#${this.anchor} div { cursor: pointer !important; }`;
      document.head.appendChild(style);
    }
    // this.chart persists across re-renders (echarts-plotting.js's renderPlot() only
    // creates it once), so the click listener must be bound once too, not on every
    // afterPlotCreated() call (which would stack duplicate handlers).
    if (this._chartClickBound) return;
    this._chartClickBound = true;
    this.chart.on("click", (params) => {
      if (params.componentType !== "series" || !Array.isArray(params.data)) return;
      this.openDrilldown(params.data[2]);
    });
  }
}

// Make EchartsSeqContentPlot globally available (referenced by bare
// window.EchartsSeqContentPlot in echarts-plotting.js's initPlot(), not a bare identifier,
// so it's robust across the minified bundle regardless of module ordering after bundling).
window.EchartsSeqContentPlot = EchartsSeqContentPlot;
