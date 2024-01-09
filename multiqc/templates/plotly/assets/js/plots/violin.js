class ViolinPlot extends Plot {
  constructor(dump) {
    super(dump);
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.active_dataset_idx].data_by_metric.length; // no samples in a dataset
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let data_by_metric = this.datasets[this.active_dataset_idx].data_by_metric;
    let header_by_metric = this.datasets[this.active_dataset_idx].header_by_metric;
    let samples = this.datasets[this.active_dataset_idx].samples;
    if (data_by_metric.length === 0) return [];

    let layout = this.layout;

    return Object.keys(header_by_metric).map((metric, idx) => {
      let params = JSON.parse(JSON.stringify(this.trace_params)); // deep copy

      let key = idx === 0 ? "" : idx + 1;
      layout["xaxis" + key] = Object.assign(
        JSON.parse(JSON.stringify(layout.xaxis)),
        header_by_metric[metric]["xaxis"],
      );
      layout["yaxis" + key] = JSON.parse(JSON.stringify(layout.yaxis));

      let sampleIds = [];
      Object.keys(data_by_metric[metric]).map((sample) => {
        sampleIds.push(samples.indexOf(sample));
      });

      return {
        type: "violin",
        x: Object.values(data_by_metric[metric]),
        name: header_by_metric[metric].title + "  ",
        text: Object.keys(data_by_metric[metric]), // sample names
        customdata: sampleIds,
        fillcolor: header_by_metric[metric].color,
        xaxis: "x" + key,
        yaxis: "y" + key,
        ...params,
        hoveron: "points",
      };
    });
  }

  afterPlotCreated() {
    let target = this.target;
    let plot = document.getElementById(target);

    let metrics = Object.keys(this.datasets[this.active_dataset_idx].header_by_metric);

    plot.on("plotly_hover", function (e) {
      let point = e.points[0];
      let pointNum = point.pointNumber;

      console.log("hovered", point, pointNum);

      let allPointsToHover = metrics.map((metric, curveNum) => {
        return { curveNumber: curveNum, pointNumber: pointNum };
      });
      let allAxisToHover = metrics.map((metric, curveNum) => {
        let key = curveNum === 0 ? "" : curveNum + 1;
        return "x" + key + "y" + key;
      });

      console.log("triggering hover of points", allPointsToHover, allAxisToHover);
      Plotly.Fx.hover(target, allPointsToHover, allAxisToHover);
    });
  }
}
