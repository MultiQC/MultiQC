// Build an rgba() color string at the given alpha from a category color.
// cat.color (multiqc/utils/mqc_colour.py::color_to_rgb_string) is always already a
// full "rgb(r,g,b)" or "rgba(r,g,b,a)" CSS string, never bare "r,g,b" digits (verified
// against a live report: FastQC's default-palette cats serialize as e.g.
// "rgba(124,181,236,1)"). Naively concatenating "rgba(" + cat.color + "," + alpha + ")"
// (the pattern default template's Plotly bar.js uses) therefore nests the color inside
// a second rgba() wrapper, which is invalid CSS; extract the r,g,b components instead
// and always emit a fresh, valid rgba().
function colorWithAlpha(color, alpha) {
  let match = color.match(/rgba?\(([^)]+)\)/);
  let [r, g, b] = (match ? match[1] : color).split(",").map((s) => s.trim());
  return `rgba(${r},${g},${b},${alpha})`;
}

// ECharts bar plot: extends the default template's BarPlot (imported just before this
// file in main-js.js) purely for its prepData()/exportData()/formatDatasetForAiPrompt()/
// activeDatasetSize(), which are engine-neutral (they only read/filter raw dataset
// fields through applyToolboxSettings()). window.BarPlot itself extends window.Plot,
// which by this point in the import order is OUR echarts Plot class (from
// echarts-plotting.js), not Plotly's.
class EchartsBarPlot extends window.BarPlot {
  // Constructs and returns ECharts bar series for the active dataset, applying the
  // toolbox (hide/rename/highlight) and pct toggle via the inherited prepData().
  buildSeries() {
    let dataset = this.datasets[this.activeDatasetIdx];

    // Multicategory (sample groups) is out of scope for Phase 0 (see BUILD_PLAN.md
    // Task 0.5); default's prepData() would silently take its Plotly-only multicategory
    // branch instead, so we crash loudly here rather than let that happen.
    if (dataset["group_labels"]) {
      throw new Error("ECharts bar plot does not support sample groups yet");
    }

    let [cats] = this.prepData();

    // "N samples hidden" banner + hiding the whole plot group when every sample is
    // hidden: plotting-shared.js's applyToolboxSettings() has this exact logic
    // (keyed off a `plotAnchor` argument), but no BarPlot caller (default template's
    // included) ever passes that argument, so it's dead code in both templates today.
    // Rather than thread plotAnchor through shared/default code (out of scope: "Do NOT
    // edit the default template"), recompute the same signal here from data already on
    // this instance. This runs on every render, so it is naturally correct and reverses
    // itself when the toolbox hide filters are cleared.
    this._updateHiddenSamplesWarning(dataset["samples"].length, this.filteredSettings.length);

    if (cats.length === 0 || this.filteredSettings.length === 0) {
      this._axisData = { axis: "yAxis", data: [] };
      return [];
    }

    // Stashed here for renderPlot() to assign to option.yAxis.data (samples are
    // toolbox-dependent, so the skeleton from Python never contains them).
    this._axisData = { axis: "yAxis", data: this.filteredSettings.map((s) => s.name) };

    let highlighted = this.filteredSettings.filter((s) => s.highlight);
    let barmode = this.echarts.datasets[this.activeDatasetIdx].layout["_mqc"]?.barmode;
    let isGroup = barmode === "group";

    // parity gap (Phase 3): the Plotly bar plot also colors the sample-name axis TICK
    // LABELS for highlighted samples (recalculateTicks() writing tickvals/ticktext HTML
    // spans, templates/default/src/js/plots/bar.js). ECharts ignores tickvals/ticktext;
    // an equivalent would need yAxis.axisLabel.formatter + a `rich` style map keyed per
    // sample. Bar dimming below is the primary highlight signal and matches Plotly;
    // tick-label coloring is deferred.
    return cats.map((cat) => {
      let series = {
        type: "bar",
        name: cat.name,
        barCategoryGap: "30%",
        data: this.filteredSettings.map((sample, sampleIdx) => {
          // Dim non-highlighted samples to 0.1 alpha when any highlight is active,
          // same rule as default plots/bar.js.
          let alpha = highlighted.length > 0 && sample.highlight === null ? 0.1 : 1;
          return {
            value: cat.data[sampleIdx],
            itemStyle: { color: colorWithAlpha(cat.color, alpha) },
          };
        }),
      };
      if (!isGroup) series.stack = "total";
      return series;
    });
  }

  // Tooltip formatter can't live in the JSON-safe serialized skeleton, so it's attached
  // here to the fully-assembled option (see Plot.applyOptionOverrides in
  // echarts-plotting.js). Without this, ECharts' default axis-tooltip prints raw float
  // values (no rounding).
  applyOptionOverrides(option) {
    if (!option.tooltip) return;
    option.tooltip.formatter = (params) => {
      let list = Array.isArray(params) ? params : [params];
      if (list.length === 0) return "";
      let sampleLabel = list[0].name;
      let rows = list.map((p) => `${p.marker}${p.seriesName}: <b>${window.formatNumber(p.value)}</b>`).join("<br/>");
      return `${sampleLabel}<br/>${rows}`;
    };
  }

  // See the comment at the buildSeries() call site above.
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
}

// Make EchartsBarPlot globally available (referenced by bare window.EchartsBarPlot in
// echarts-plotting.js's initPlot(), not a bare identifier, so it's robust across the
// minified bundle regardless of module ordering after bundling).
window.EchartsBarPlot = EchartsBarPlot;
