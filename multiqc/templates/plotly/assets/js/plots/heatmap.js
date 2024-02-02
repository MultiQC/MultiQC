class HeatmapPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.xcats_samples = dump["xcats_samples"];
    this.ycats_samples = dump["ycats_samples"];
    this.square = dump["square"];
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    let rows = this.datasets[this.activeDatasetIdx]["rows"];
    if (rows.length === 0) return 0; // no rows
    return rows[0].length; // no columns in a row
  }

  prepData() {
    // Prepare data to either build Plotly traces or export as a file
    let dataset = this.datasets[this.activeDatasetIdx];
    let rows = dataset["rows"];
    let xcats = dataset["xcats"];
    let ycats = dataset["ycats"];

    if (this.xcats_samples) {
      let xcatsSettings = applyToolboxSettings(xcats);
      if (xcatsSettings == null) return; // All series are hidden, do not render the graph.

      // If sample is highlighted, make all values None except for the highlighted sample
      // But only if it's not done for the Y
      if (!this.ycats_samples) {
        if (xcatsSettings.filter((s) => s.highlight).length > 0) {
          rows = rows.map((row) => {
            return row.map((val, x) => (xcatsSettings[x].highlight ? val : null));
          });
        }
      }
      // Rename and filter samples:
      xcats = xcatsSettings.filter((s) => !s.hidden).map((s) => s.name);
      rows = rows.map((row) => row.filter((val, i) => !xcatsSettings[i].hidden));
    }

    if (this.ycats_samples) {
      let ycatsSettings = applyToolboxSettings(ycats);
      if (ycatsSettings == null) return; // All series are hidden, do not render the graph.

      // If sample is highlighted, make all values None except for the highlighted sample
      if (ycatsSettings.filter((s) => s.highlight).length > 0) {
        rows = rows.map((row, y) => {
          return ycatsSettings[y].highlight ? row : row.map(() => null);
        });
      }
      // Rename and filter samples:
      ycats = ycatsSettings.filter((s) => !s.hidden).map((s) => s.name);
      rows = rows.filter((row, i) => !ycatsSettings[i].hidden);
    }

    return [rows, xcats, ycats];
  }

  buildTraces(layout) {
    let [rows, xcats, ycats] = this.prepData();

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
    let dataset = this.datasets[this.activeDatasetIdx];
    let xcats = dataset["xcats"];
    let ycats = dataset["ycats"];

    let pxPerElem = newHeight / ycats.length;
    let newWidth = pxPerElem * xcats.length;
    newWidth = this.square ? newWidth : null;

    if (newHeight < this.layout.height) {
      // We're shrinking the plot, so we need to allow plotly to skip ticks
      dataset.layout.xaxis.nticks = null;
      dataset.layout.yaxis.nticks = null;
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
