// ECharts heatmap plot: STANDALONE class extending window.Plot directly. Unlike
// bar/line/scatter, this file does NOT import the default template's plots/heatmap.js:
// that file's top-level `$(function () {...})` handler (bottom of the file) calls
// `Plotly.restyle(anchor, {zmin: ...})` on the zmin/zmax range-slider inputs and toggles
// clustering, and would run (and crash on the undefined `Plotly` global) as soon as this
// bundle loads. Instead, `prepData()`/`exportData()`/`formatDatasetForAiPrompt()` are
// copied field-for-field from that file, and the two bottom handlers are re-implemented
// below for ECharts (visualMap min/max instead of Plotly.restyle; renderPlot is already
// engine-neutral, just copied here since it lives in the file we don't import).
class EchartsHeatmapPlot extends window.Plot {
  constructor(dump) {
    super(dump);
    this.xCatsAreSamples = dump["xcats_samples"];
    this.yCatsAreSamples = dump["ycats_samples"];
    this.filtXCatsSettings = [];
    this.filtYCatsSettings = [];
    this.clusterSwitchClusteredActive = dump["cluster_switch_clustered_active"];
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    let rows = this.datasets[this.activeDatasetIdx]["rows"];
    if (rows.length === 0) return 0; // no rows
    return rows[0].length; // no columns in a row
  }

  // Copied from templates/default/src/js/plots/heatmap.js: resolves the clustered vs.
  // unclustered rows/xcats/ycats for the current cluster-switch state, then applies
  // toolbox hide/rename/highlight to whichever axis carries sample names.
  prepData(dataset) {
    dataset = dataset ?? this.datasets[this.activeDatasetIdx];
    let rows =
      this.clusterSwitchClusteredActive && dataset["rows_clustered"] ? dataset["rows_clustered"] : dataset["rows"];
    let xcats =
      this.clusterSwitchClusteredActive && dataset["xcats_clustered"] ? dataset["xcats_clustered"] : dataset["xcats"];
    let ycats =
      this.clusterSwitchClusteredActive && dataset["ycats_clustered"] ? dataset["ycats_clustered"] : dataset["ycats"];

    if (this.xCatsAreSamples) {
      let xcatsSettings = applyToolboxSettings(xcats);
      if (xcatsSettings === null) return [[], [], []]; // all series hidden

      rows = rows.map((row) => row.filter((val, i) => !xcatsSettings[i].hidden));
      this.filtXCatsSettings = xcatsSettings.filter((s) => !s.hidden);
      xcats = this.filtXCatsSettings.map((s) => s.name);
    } else {
      this.filtXCatsSettings = [];
    }

    if (this.yCatsAreSamples) {
      let yCatsSettings = applyToolboxSettings(ycats);
      if (yCatsSettings === null) return [[], [], []]; // all series hidden

      rows = rows.filter((row, i) => !yCatsSettings[i].hidden);
      this.filtYCatsSettings = yCatsSettings.filter((s) => !s.hidden);
      ycats = this.filtYCatsSettings.map((s) => s.name);
    } else {
      this.filtYCatsSettings = [];
    }

    return [rows, xcats, ycats];
  }

  plotAiHeader() {
    let result = super.plotAiHeader();
    if (this.pconfig.xlab) result += `X axis: ${this.pconfig.xlab}\n`;
    if (this.pconfig.ylab) result += `Y axis: ${this.pconfig.ylab}\n`;
    if (this.pconfig.zlab) result += `Z axis: ${this.pconfig.zlab}\n`;
    return result;
  }

  // Copied verbatim (same field access) from default heatmap.js.
  formatDatasetForAiPrompt(dataset) {
    let prompt = "";
    let [rows, xcats, ycats] = this.prepData(dataset);

    if (xcats.length === 0 || ycats.length === 0) {
      prompt +=
        "All samples are hidden by user, so no data to analyse. Please inform user to use the toolbox to unhide samples.\n";
      return prompt;
    }

    if (xcats) {
      if (ycats) {
        prompt = "|";
        if (this.yCatsAreSamples) prompt += "Sample";
      }
      if (this.xCatsAreSamples) {
        const xPseudonyms = this.filtXCatsSettings.map((s) => s.pseudonym ?? s.name);
        prompt += "|" + xPseudonyms.join("|") + "|\n";
      } else {
        prompt += "|" + xcats.join("|") + "|\n";
      }
      if (ycats) prompt += "|---";
      prompt += "|" + xcats.map(() => "---").join("|") + "|\n";
    }
    for (let i = 0; i < rows.length; i++) {
      if (ycats) {
        if (this.yCatsAreSamples) {
          const yPseudonyms = this.filtYCatsSettings.map((s) => s.pseudonym ?? s.name);
          prompt += "|" + yPseudonyms[i];
        } else {
          prompt += "|" + ycats[i];
        }
      }
      prompt +=
        "|" +
        rows[i].map((x) => (!Number.isFinite(x) ? "" : Number.isInteger(x) ? x : parseFloat(x.toFixed(2)))).join("|") +
        "|\n";
    }
    return prompt;
  }

  // Copied verbatim from default heatmap.js.
  exportData(format) {
    let [rows, xcats, ycats] = this.prepData();

    let sep = format === "tsv" ? "\t" : ",";

    let csv = [".", ...xcats].join(sep) + "\n";
    for (let i = 0; i < rows.length; i++) {
      csv += [ycats[i], ...rows[i]].join(sep) + "\n";
    }
    return csv;
  }

  // ECharts analog of default heatmap.js's buildTraces(): builds cell data
  // `[xi, yi, val]` from the toolbox-filtered rows/xcats/ycats. Cell labels are shown
  // when the Python side enabled them (`dataset.trace_params.texttemplate`, set by
  // `HeatmapPlot.create` when there are few enough cells or `pconfig.display_values`),
  // formatted to `pconfig.tt_decimals` decimal places to match the Plotly texttemplate.
  buildSeries() {
    let [rows, xcats, ycats] = this.prepData();
    if (rows.length === 0 || xcats.length === 0 || ycats.length === 0) return [];

    // Stashed for buildCurrentOption() (xAxis, generic single-axis path) and
    // applyOptionOverrides() below (yAxis, since a heatmap needs both).
    this._axisData = { axis: "xAxis", data: xcats };
    this._yAxisData = ycats;

    let dataset = this.datasets[this.activeDatasetIdx];
    let showLabels = !!(dataset["trace_params"] && dataset["trace_params"]["texttemplate"]);
    let decimals = this.pconfig.tt_decimals ?? 2;

    let data = [];
    for (let yi = 0; yi < rows.length; yi++) {
      for (let xi = 0; xi < xcats.length; xi++) {
        let val = rows[yi][xi];
        if (val === null || val === undefined) continue;
        if (showLabels && typeof val === "number") {
          data.push({ value: [xi, yi, val], label: { show: true, formatter: val.toFixed(decimals) } });
        } else {
          data.push([xi, yi, val]);
        }
      }
    }

    return [{ type: "heatmap", data: data }];
  }

  // buildCurrentOption() (echarts-plotting.js) only patches ONE axis generically from
  // `plot._axisData` (bar/line only ever need at most one toolbox-dependent category
  // axis); a heatmap needs both, so the y axis and the tooltip formatter (which can
  // never live in the JSON-safe serialized skeleton) are applied here instead.
  applyOptionOverrides(option) {
    if (option.yAxis && this._yAxisData) option.yAxis.data = this._yAxisData;

    if (option.tooltip) {
      let xlab = this.pconfig.xlab || "x";
      let ylab = this.pconfig.ylab || "y";
      let zlab = this.pconfig.zlab || "z";
      let xcats = this._axisData ? this._axisData.data : [];
      let ycats = this._yAxisData ?? [];
      option.tooltip.formatter = (params) => {
        let [xi, yi, val] = params.value;
        return `${xlab}: ${xcats[xi]}<br/>${ylab}: ${ycats[yi]}<br/>${zlab}: ${val}`;
      };
    }
  }

  resize(newHeight) {
    let dataset = this.datasets[this.activeDatasetIdx];
    let xcats = dataset["xcats"];
    let ycats = dataset["ycats"];
    let newWidth = null;
    if (this.square && ycats.length > 0) {
      let pxPerElem = (newHeight - 200) / ycats.length;
      newWidth = pxPerElem * xcats.length + 250;
    }
    super.resize(newHeight, newWidth);
  }
}

// Make EchartsHeatmapPlot globally available (referenced by bare window.EchartsHeatmapPlot
// in echarts-plotting.js's initPlot(), not a bare identifier, so it's robust across the
// minified bundle regardless of module ordering after bundling).
window.EchartsHeatmapPlot = EchartsHeatmapPlot;

$(function () {
  // zmin/zmax range sliders: re-implemented for ECharts as a live visualMap update
  // instead of default heatmap.js's `Plotly.restyle(anchor, {zmin: ...})`.
  $(".mqc_hcplot_range_sliders input").on("keyup change input", function () {
    let anchor = $(this).data("anchor");
    let minmax = $(this).data("minmax"); // "min" or "max"
    let plot = mqc_plots[anchor];
    let value = parseFloat($(this).val());
    if (plot && plot.chart && !Number.isNaN(value)) {
      plot.chart.setOption({ visualMap: { [minmax]: value } });
    }
    $("#" + anchor + "_range_slider_" + minmax + ", #" + anchor + "_range_slider_" + minmax + "_txt").val(
      $(this).val(),
    );
  });

  // Cluster toggle: engine-neutral (set a flag, then re-render from scratch), copied
  // here (rather than imported) only because it lives in the default heatmap.js file we
  // deliberately don't import (see the file header comment).
  $('button[data-action="unclustered"], button[data-action="clustered"]').on("click", function (e) {
    e.preventDefault();
    let $btn = $(this);
    let plotAnchor = $(this).data("plot-anchor");
    let plot = mqc_plots[plotAnchor];
    if (!plot) return;

    // Toggle buttons
    $btn.toggleClass("active").siblings().toggleClass("active");

    // Update plot
    plot.clusterSwitchClusteredActive = $btn.data("action") === "clustered";
    renderPlot(plotAnchor); // re-render
  });
});
