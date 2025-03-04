class HeatmapPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.xCatsAreSamples = dump["xcats_samples"];
    this.yCatsAreSamples = dump["ycats_samples"];
    this.square = dump["square"];
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

  prepData(dataset) {
    // Prepare data to either build Plotly traces or export as a file
    dataset = dataset ?? this.datasets[this.activeDatasetIdx];
    let rows =
      this.clusterSwitchClusteredActive && dataset["rows_clustered"] ? dataset["rows_clustered"] : dataset["rows"];
    let xcats =
      this.clusterSwitchClusteredActive && dataset["xcats_clustered"] ? dataset["xcats_clustered"] : dataset["xcats"];
    let ycats =
      this.clusterSwitchClusteredActive && dataset["ycats_clustered"] ? dataset["ycats_clustered"] : dataset["ycats"];

    if (this.xCatsAreSamples) {
      let xcatsSettings = applyToolboxSettings(xcats);
      if (xcatsSettings === null) return; // All series are hidden, do not render the graph.

      rows = rows.map((row) => row.filter((val, i) => !xcatsSettings[i].hidden));
      this.filtXCatsSettings = xcatsSettings.filter((s) => !s.hidden);
      xcats = this.filtXCatsSettings.map((s) => s.name);
    }

    if (this.yCatsAreSamples) {
      let yCatsSettings = applyToolboxSettings(ycats);
      if (yCatsSettings === null) return; // All series are hidden, do not render the graph.

      rows = rows.filter((row, i) => !yCatsSettings[i].hidden);
      this.filtYCatsSettings = yCatsSettings.filter((s) => !s.hidden);
      ycats = this.filtYCatsSettings.map((s) => s.name);
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

  formatDatasetForAiPrompt(dataset) {
    let prompt = "";
    let [rows, xcats, ycats] = this.prepData(dataset);

    // Check if all samples are hidden
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

  buildTraces() {
    let [rows, xcats, ycats] = this.prepData();
    if (rows.length === 0 || xcats.length === 0 || ycats.length === 0) return [];

    if (this.filtYCatsSettings.length > 0) {
      const maxYTicks = (this.layout.height - 200) / 12;
      this.recalculateTicks(this.filtYCatsSettings, this.layout.yaxis, maxYTicks);
    }
    if (this.filtXCatsSettings.length > 0) {
      const maxXTicks = (this.layout.width - 250) / 18;
      this.recalculateTicks(this.filtXCatsSettings, this.layout.xaxis, maxXTicks);
    }

    let dataset = this.datasets[this.activeDatasetIdx];
    let params = JSON.parse(JSON.stringify(dataset["trace_params"])); // deep copy

    let heatmap = {
      type: "heatmap",
      z: rows,
      x: xcats,
      y: ycats,
      ...params,
    };
    return [heatmap];
  }

  exportData(format) {
    let [rows, xcats, ycats] = this.prepData();

    let sep = format === "tsv" ? "\t" : ",";

    let csv = [".", ...xcats].join(sep) + "\n";
    for (let i = 0; i < rows.length; i++) {
      csv += [ycats[i], ...rows[i]].join(sep) + "\n";
    }
    return csv;
  }

  resize(newHeight) {
    this.layout.height = newHeight;

    const maxYTicks = (this.layout.height - 200) / 12;
    this.recalculateTicks(this.filtYCatsSettings, this.layout.yaxis, maxYTicks);

    let dataset = this.datasets[this.activeDatasetIdx];
    let xcats = dataset["xcats"];
    let ycats = dataset["ycats"];
    let pxPerElem = (newHeight - 200) / ycats.length;
    let newWidth = null;
    if (this.square) newWidth = pxPerElem * xcats.length + 250;

    super.resize(newHeight, newWidth);
  }
}

$(function () {
  // Listeners for range slider
  $(".mqc_hcplot_range_sliders input").on("keyup change input", function () {
    let anchor = $(this).data("anchor");
    let minmax = $(this).data("minmax");
    if (minmax === "min") Plotly.restyle(anchor, { zmin: $(this).val() });
    if (minmax === "max") Plotly.restyle(anchor, { zmax: $(this).val() });
    $("#" + anchor + "_range_slider_" + minmax + ", #" + anchor + "_range_slider_" + minmax + "_txt").val(
      $(this).val(),
    );
  });

  // Listener for clustering toggle
  $('button[data-action="unclustered"], button[data-action="clustered"]').on("click", function (e) {
    e.preventDefault();
    let $btn = $(this);
    let plotAnchor = $(this).data("plot-anchor");
    let plot = mqc_plots[plotAnchor];

    // Toggle buttons
    $btn.toggleClass("active").siblings().toggleClass("active");

    // Update plot
    plot.clusterSwitchClusteredActive = $btn.data("action") === "clustered";
    renderPlot(plotAnchor); // re-render
  });
});
