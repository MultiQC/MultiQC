class ViolinPlot extends Plot {
  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.active_dataset_idx].data_by_metric.length; // no samples in a dataset
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let data_by_metric = this.datasets[this.active_dataset_idx].data_by_metric;
    let header_by_metric = this.datasets[this.active_dataset_idx].header_by_metric;
    if (data_by_metric.length === 0) return [];

    this.layout.grid = {
      rows: Object.keys(data_by_metric).length,
      columns: 1,
      pattern: "independent",
      roworder: "top to bottom",
      ygap: 0.3,
    };

    let layout = this.layout;
    layout.margin.l = 160;
    layout.margin.pad = 20;

    return Object.keys(data_by_metric).map((metric, idx) => {
      let params = JSON.parse(JSON.stringify(this.trace_params)); // deep copy

      let akey = "xaxis" + (idx === 0 ? "" : idx + 1);
      let axis = layout[akey] ?? {};
      axis.rangemode = header_by_metric[metric]["tozero"] ? "tozero" : "normal";
      axis.range = [header_by_metric[metric]["min"], header_by_metric[metric]["max"]];
      layout[akey] = axis;

      return {
        type: "violin",
        x: Object.values(data_by_metric[metric]),
        name: header_by_metric[metric].title,
        text: Object.keys(data_by_metric[metric]), // sample names
        xaxis: "x" + idx,
        yaxis: "y" + idx,
        range: [0, 100],
        ...params,
      };
    });
  }
}
