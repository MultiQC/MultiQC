class HeatmapPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.xcats_samples = dump.xcats_samples;
    this.ycats_samples = dump.ycats_samples;
    this.square = dump.square;
    this.heatmap_config = dump.heatmap_config;
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    let rows = this.datasets[this.active_dataset_idx].rows;
    if (rows.length === 0) return 0; // no rows
    return rows[0].length; // no columns in a row
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let rows = this.datasets[this.active_dataset_idx].rows;
    let xcats = this.datasets[this.active_dataset_idx].xcats;
    let ycats = this.datasets[this.active_dataset_idx].ycats;

    if (this.xcats_samples) {
      let xcatsSettings = applyToolboxSettings(xcats);
      if (xcatsSettings == null) return; // All series are hidden, do not render the graph.
      // Rename and filter samples:
      xcats = xcatsSettings.filter((s) => !s.hidden).map((s) => s.name);
      rows = rows.map((row) => row.filter((val, i) => !xcatsSettings[i].hidden));
    }

    if (this.ycats_samples) {
      let ycatsSettings = applyToolboxSettings(ycats);
      if (ycatsSettings == null) return; // All series are hidden, do not render the graph.
      // Rename and filter samples:
      ycats = ycatsSettings.filter((s) => !s.hidden).map((s) => s.name);
      rows = rows.filter((row, i) => !ycatsSettings[i].hidden);
    }

    let heatmap = {
      type: "heatmap",
      z: rows,
      x: xcats,
      y: ycats,
    };
    heatmap = Object.assign(heatmap, this.heatmap_config);
    return [heatmap];
  }

  resize(newHeight) {
    let xcats = this.datasets[this.active_dataset_idx].xcats;
    let ycats = this.datasets[this.active_dataset_idx].ycats;

    let pxPerElem = newHeight / ycats.length;
    let newWidth = pxPerElem * xcats.length;
    newWidth = this.square ? newWidth : null;

    if (newHeight < this.layout.height) {
      // We're shrinking the plot, so we need to allow plotly to skip ticks
      this.layout.xaxis.nticks = null;
      this.layout.yaxis.nticks = null;
    }

    Plotly.relayout(this.target, { height: newHeight, width: newWidth });
  }
}

$(function () {
  // Listeners for range slider
  $(".mqc_hcplot_range_sliders input").on("keyup change input", function () {
    let target = $(this).data("target");
    let minmax = $(this).data("minmax");
    if (minmax === "min") Plotly.restyle(target, { zmin: $(this).val() });
    if (minmax === "max") Plotly.restyle(target, { zmax: $(this).val() });
    $("#" + target + "_range_slider_" + minmax + ", #" + target + "_range_slider_" + minmax + "_txt").val(
      $(this).val(),
    );
  });
});
