class ScatterPlot extends Plot {
  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx].points; // no data points in a dataset
  }

  prepData() {
    // Prepare data to either build Plotly traces or export as a file
    let dataset = this.datasets[this.activeDatasetIdx];

    let points = dataset.points;

    let samples = points.map((point) => point.name);
    let sampleSettings = applyToolboxSettings(samples);

    points = points.map((point, idx) => {
      point.name = sampleSettings[idx].name ?? point.name;
      point.highlight = sampleSettings[idx].highlight;
      if (!sampleSettings[idx].hidden) return point;
    });

    return [samples, points];
  }

  buildTraces() {
    let dataset = this.datasets[this.activeDatasetIdx];

    let [samples, points] = this.prepData();
    if (points.length === 0 || samples.length === 0) return [];

    // Reorder points so highlighted points are on top
    let highlighted = points.filter((p) => p.highlight);
    let nonHighlighted = points.filter((p) => !p.highlight);
    points = nonHighlighted.concat(highlighted);

    return points.map((point) => {
      let params = JSON.parse(JSON.stringify(dataset["trace_params"])); // deep copy
      params.marker.size = point["marker_size"] ?? params.marker.size;
      params.marker.line = {
        width: point["marker_line_width"] ?? params.marker.line.width,
      };
      params.marker.opacity = point["opacity"] ?? params.marker.opacity;
      params.marker.color = point["color"] ?? params.marker.color;
      if (highlighted.length > 0) params.marker.color = point.highlight ?? "#cccccc";

      return {
        type: "scatter",
        x: [point.x],
        y: [point.y],
        name: point.name,
        text: [point.annotation ?? point.name],
        ...params,
      };
    });
  }

  exportData(format) {
    let [samples, points] = this.prepData();

    let sep = format === "tsv" ? "\t" : ",";

    let csv = ["Name", "X", "Y"].join(sep) + "\n";
    for (let i = 0; i < samples.length; i++) {
      let point = points[i];
      csv += [samples[i], +point.x, point.y].join(sep) + "\n";
    }
    return csv;
  }
}
