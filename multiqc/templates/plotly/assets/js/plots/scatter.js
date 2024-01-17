class ScatterPlot extends Plot {
  constructor(data) {
    super(data);
    this.categories = data.categories;
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx].points; // no data points in a dataset
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let points = this.datasets[this.activeDatasetIdx].points;
    if (points.length === 0) return [];

    let samples = points.map((point) => point.name);
    let sampleSettings = applyToolboxSettings(samples);
    if (sampleSettings == null) return; // All series are hidden, do not render the graph
    points = points.map((point, idx) => {
      point.name = sampleSettings[idx].name ?? point.name;
      point.highlight = sampleSettings[idx].highlight;
      if (!sampleSettings[idx].hidden) return point;
    });

    // Reorder points so highlighted points are on top
    let highlighted = points.filter((p) => p.highlight);
    let nonHighlighted = points.filter((p) => !p.highlight);
    points = nonHighlighted.concat(highlighted);

    return points.map((point) => {
      let x = point.x;
      if (this.categories && Number.isInteger(x) && x < this.categories.length) x = this.categories[x];

      let params = JSON.parse(JSON.stringify(this.traceParams)); // deep copy
      params.marker.size = point["marker_size"] ?? params.marker.size;
      params.marker.line = {
        width: point["marker_line_width"] ?? params.marker.line.width,
      };
      params.marker.opacity = point["opacity"] ?? params.marker.opacity;
      params.marker.color = point["color"] ?? params.marker.color;
      if (highlighted.length > 0) params.marker.color = point.highlight ?? "#cccccc";

      return {
        type: "scatter",
        x: [x],
        y: [point.y],
        name: point.name,
        text: [point.annotation ?? point.name],
        ...params,
      };
    });
  }
}
