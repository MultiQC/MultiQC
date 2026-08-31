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
// `[start, end, row, r, g, b]` data item to the same rect geometry and fill color; a bin
// covers columns `start..end` inclusive (1-based), so its rect spans the data-space interval
// `[start, end + 1)` on the value xAxis, and `row` (a category-axis index) is used as-is.
function seqContentRenderItem(params, api) {
  var start = api.value(0);
  var end = api.value(1);
  var row = api.value(2);
  var r = api.value(3);
  var g = api.value(4);
  var b = api.value(5);
  var p0 = api.coord([start, row]);
  var p1 = api.coord([end + 1, row]);
  var height = api.size([0, 1])[1];
  return {
    type: "rect",
    shape: { x: p0[0], y: p0[1] - height / 2, width: p1[0] - p0[0], height: height },
    style: { fill: "rgb(" + r + "," + g + "," + b + ")" },
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

    let data = [];
    for (let rowIdx = 0; rowIdx < rows.length; rowIdx++) {
      for (let bin of rows[rowIdx]) {
        let [r, g, b] = seqContentBinRgb(bin);
        data.push([bin.start, bin.end, rowIdx, r, g, b]);
      }
    }

    return [
      {
        type: "custom",
        coordinateSystem: "cartesian2d",
        renderItem: seqContentRenderItem,
        data: data,
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
    option.tooltip.formatter = (params) => {
      let [start, end, rowIdx] = params.value;
      let bin = (rows[rowIdx] ?? []).find((b) => b.start === start && b.end === end);
      if (!bin) return "";
      let sampleName = axisData[rowIdx] ?? "";
      let posLabel = start === end ? `${start} bp` : `${start}-${end} bp`;
      return (
        `<b>${sampleName}</b><br/>${posLabel}<br/>` +
        `%T: ${bin.t}%<br/>%C: ${bin.c}%<br/>%A: ${bin.a}%<br/>%G: ${bin.g}%`
      );
    };
  }

  // No click/drilldown here (Phase 3); afterPlotCreated() left as the inherited no-op
  // from echarts-plotting.js's Plot base class, a clean seam for T3.1 to fill in.
}

// Make EchartsSeqContentPlot globally available (referenced by bare
// window.EchartsSeqContentPlot in echarts-plotting.js's initPlot(), not a bare identifier,
// so it's robust across the minified bundle regardless of module ordering after bundling).
window.EchartsSeqContentPlot = EchartsSeqContentPlot;
