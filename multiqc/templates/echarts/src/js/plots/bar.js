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

    if (dataset["group_labels"]) {
      return this.buildGroupedSeries(dataset);
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

  // Grouped ("sample groups" / multicategory) bars, e.g. riboWaltz's per-region frame
  // plots. `dataset.group_labels`/`dataset.offset_groups` are the same fields the
  // default template's BarPlot.prepData()/buildTraces() multicategory branch reads
  // (templates/default/src/js/plots/bar.js); the inherited prepData() already does the
  // group-aware toolbox work (hide/rename/highlight BY GROUP, not by sample) and sets
  // this.filteredSettings/originalGroupLabels/filteredGroupLabels/offsetGroups/
  // groupSettingsMap accordingly, filtered to visible groups.
  //
  // Plotly draws this with a shared `offsetgroup` per sample on its native
  // multicategory axis. ECharts has no multicategory axis, so the pragmatic port is:
  // unique group labels become the yAxis categories, and each sample becomes one
  // ECharts `stack` id positioned at its group's slot. Categories within one sample's
  // stack sum together (mirrors Plotly's overlapping same-offsetgroup traces); ECharts
  // automatically places different stack ids side by side at a shared category, which
  // is what makes the groups visually distinct.
  buildGroupedSeries(dataset) {
    let [cats] = this.prepData();

    let totalGroups = new Set(dataset["group_labels"]).size;
    let visibleGroups = new Set(this.filteredGroupLabels).size;
    this._updateHiddenSamplesWarning(totalGroups, visibleGroups, "groups");

    if (cats.length === 0 || this.filteredSettings.length === 0) {
      this._axisData = { axis: "yAxis", data: [] };
      return [];
    }

    let groupAxis = [];
    let groupIndex = {};
    this.filteredGroupLabels.forEach((label) => {
      if (!(label in groupIndex)) {
        groupIndex[label] = groupAxis.length;
        groupAxis.push(label);
      }
    });
    this._axisData = { axis: "yAxis", data: groupAxis };

    // Rows sharing the same (row-level) sample name belong to the same stack, e.g.
    // riboWaltz's "sample1_28nt"/"sample1_29nt" rows both display as "sample1" and
    // must stack together to form one compound bar per length group.
    let rowsBySample = {};
    this.filteredSettings.forEach((row, rowIdx) => {
      if (!rowsBySample[row.name]) rowsBySample[row.name] = [];
      rowsBySample[row.name].push(rowIdx);
    });

    let anyGroupHighlight = Object.values(this.groupSettingsMap || {}).some((g) => g.highlight);

    let series = [];
    Object.keys(rowsBySample).forEach((sampleName) => {
      let rowIndices = rowsBySample[sampleName];
      let stackId = (this.offsetGroups && this.offsetGroups[sampleName]) || sampleName;

      cats.forEach((cat) => {
        let data = new Array(groupAxis.length).fill(null);
        rowIndices.forEach((rowIdx) => {
          let originalGroupLabel = this.originalGroupLabels[rowIdx];
          let groupHighlight = this.groupSettingsMap?.[originalGroupLabel]?.highlight;
          let alpha = anyGroupHighlight && !groupHighlight ? 0.1 : 1;
          data[groupIndex[this.filteredGroupLabels[rowIdx]]] = {
            value: cat.data[rowIdx],
            name: sampleName, // read by the tooltip formatter to disambiguate the stack
            itemStyle: { color: colorWithAlpha(cat.color, alpha) },
          };
        });
        series.push({
          type: "bar",
          name: cat.name,
          stack: stackId,
          barCategoryGap: "30%",
          data,
        });
      });
    });
    return series;
  }

  // Tooltip formatter can't live in the JSON-safe serialized skeleton, so it's attached
  // here to the fully-assembled option (see Plot.applyOptionOverrides in
  // echarts-plotting.js). Without this, ECharts' default axis-tooltip prints raw float
  // values (no rounding).
  applyOptionOverrides(option) {
    if (!option.tooltip) return;
    option.tooltip.formatter = (params) => {
      let list = Array.isArray(params) ? params : [params];
      // Grouped bars (buildGroupedSeries) fill unused axis slots with null so a
      // sample's stack only shows up at the group(s) it belongs to; drop those from
      // the tooltip instead of printing a blank/NaN row for every other stack.
      list = list.filter((p) => p.value !== null && p.value !== undefined);
      if (list.length === 0) return "";
      let sampleLabel = list[0].name;
      let rows = list
        .map((p) => {
          // Grouped bars stash the row-level sample name on the data item (buildGroupedSeries)
          // to disambiguate which stack a row belongs to, since several samples' series
          // share the same seriesName (the category name).
          let stackLabel = p.data && p.data.name ? `${p.data.name} - ` : "";
          return `${p.marker}${stackLabel}${p.seriesName}: <b>${window.formatNumber(p.value)}</b>`;
        })
        .join("<br/>");
      return `${sampleLabel}<br/>${rows}`;
    };
  }

  // See the comment at the buildSeries() call site above.
  _updateHiddenSamplesWarning(total, visible, label = "samples") {
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
        ⚠ <strong>Warning:</strong> ${nHidden} ${label} hidden.
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
