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
    // axisLabel (item B): colors the sample-name axis tick labels to match Plotly's
    // recalculateTicks() (ticktext HTML spans); see Plot.sampleAxisLabel() in
    // echarts-plotting.js. maxTicks mirrors the Plotly bar plot's own formula
    // (templates/plotly/src/js/plots/bar.js), which derives it from the plot height.
    let maxTicks = (this.layout.height - 140) / 12;
    this._axisData = {
      axis: "yAxis",
      data: this.filteredSettings.map((s) => s.name),
      axisLabel: this.sampleAxisLabel(this.filteredSettings, maxTicks),
    };

    let barmode = this.echarts.datasets[this.activeDatasetIdx].layout["_mqc"]?.barmode;
    let isGroup = barmode === "group";

    // Highlight, item A/C: dim every non-highlighted sample once ANY highlight is
    // active anywhere (global flag, not just matches on this plot -- matches Plotly's
    // intended cross-plot behavior), and outline a highlighted sample's own bar segments
    // in its group color (Highcharts-era behavior, item C).
    let series = cats.map((cat) => {
      let s = {
        type: "bar",
        name: cat.name,
        barCategoryGap: "30%",
        data: this.filteredSettings.map((sample, sampleIdx) => {
          let alpha = window.mqc_highlight_f_texts.length > 0 && sample.highlight === null ? 0.1 : 1;
          let itemStyle = { color: colorWithAlpha(cat.color, alpha) };
          if (sample.highlight) {
            itemStyle.borderColor = sample.highlight;
            itemStyle.borderWidth = 2;
          }
          return {
            value: cat.data[sampleIdx],
            itemStyle,
          };
        }),
      };
      if (!isGroup) s.stack = "total";
      return s;
    });

    // Bar's category (sample) yAxis makes y_bands/y_lines meaningless, so the Python
    // side (multiqc/plots/echarts/bar.py) only stashes x_bands/x_lines here.
    series.push(...this.bandsLinesSeries());
    return series;
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
    // Highlight tick labels (item B): a grouped bar has no per-sample axis tick (samples
    // are side-by-side stacks); its highlightable axis IS the group-label axis, and the
    // toolbox highlights BY GROUP (groupSettingsMap). So color the group tick labels the
    // same way the simple path colors sample ticks, reusing Plot.sampleAxisLabel(). Build
    // one settings entry per axis label ({name, highlight}) from groupSettingsMap, keyed
    // by display name so it lines up with groupAxis (which carries the renamed labels).
    let groupByName = {};
    Object.values(this.groupSettingsMap || {}).forEach((g) => {
      groupByName[g.name] = g;
    });
    let groupSettings = groupAxis.map((label) => groupByName[label] || { name: label, highlight: null });
    let maxTicks = (this.layout.height - 140) / 12;
    this._axisData = { axis: "yAxis", data: groupAxis, axisLabel: this.sampleAxisLabel(groupSettings, maxTicks) };

    // Rows sharing the same (row-level) sample name belong to the same stack, e.g.
    // riboWaltz's "sample1_28nt"/"sample1_29nt" rows both display as "sample1" and
    // must stack together to form one compound bar per length group.
    let rowsBySample = {};
    this.filteredSettings.forEach((row, rowIdx) => {
      if (!rowsBySample[row.name]) rowsBySample[row.name] = [];
      rowsBySample[row.name].push(rowIdx);
    });

    let series = [];
    Object.keys(rowsBySample).forEach((sampleName) => {
      let rowIndices = rowsBySample[sampleName];
      let stackId = (this.offsetGroups && this.offsetGroups[sampleName]) || sampleName;

      cats.forEach((cat) => {
        let data = new Array(groupAxis.length).fill(null);
        rowIndices.forEach((rowIdx) => {
          let originalGroupLabel = this.originalGroupLabels[rowIdx];
          let groupHighlight = this.groupSettingsMap?.[originalGroupLabel]?.highlight;
          // Item A: dim gate keyed off the GLOBAL highlight flag, not the local
          // (this plot's) anyGroupHighlight -- once any highlight is active anywhere,
          // every non-matching group dims here too, matching Plotly's intended
          // cross-plot behavior.
          let alpha = window.mqc_highlight_f_texts.length > 0 && !groupHighlight ? 0.1 : 1;
          let itemStyle = { color: colorWithAlpha(cat.color, alpha) };
          if (groupHighlight) {
            itemStyle.borderColor = groupHighlight;
            itemStyle.borderWidth = 2;
          }
          data[groupIndex[this.filteredGroupLabels[rowIdx]]] = {
            value: cat.data[rowIdx],
            name: sampleName, // read by the tooltip formatter to disambiguate the stack
            itemStyle,
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
    series.push(...this.bandsLinesSeries());
    return series;
  }

  // Tooltip formatter can't live in the JSON-safe serialized skeleton, so it's attached
  // here to the fully-assembled option (see Plot.applyOptionOverrides in
  // echarts-plotting.js). Without this, ECharts' default axis-tooltip prints raw float
  // values (no rounding).
  applyOptionOverrides(option) {
    if (!option.tooltip) return;
    // Bar's value axis is xAxis (bars are horizontal, samples on the category yAxis),
    // matching the Plotly reference (templates/plotly/src/js/plots/bar.js's
    // `this.layout.xaxis?.hoverformat`). NOTE: unlike `ticksuffix`, `hoverformat` is a
    // Plotly-only key the neutral `AxisIR` deliberately drops (see layout.py's
    // `from_dataset_layout` docstring), so it never reaches the plot-level `this.layout`
    // for the echarts engine; it only survives on the raw per-dataset layout dict
    // (`BaseDataset.layout`, set by e.g. bargraph.py's `_dataset_layout`), which is what
    // `this.datasets[...]` exposes here.
    let hoverformat = this.datasets[this.activeDatasetIdx]?.layout?.xaxis?.hoverformat;
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
          return `${p.marker}${stackLabel}${p.seriesName}: <b>${window.formatWithHoverformat(p.value, hoverformat)}</b>`;
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
