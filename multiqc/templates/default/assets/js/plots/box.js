class BoxPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.filteredSettings = [];
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx]["samples"].length;
  }

  prepData() {
    let data = this.datasets[this.activeDatasetIdx]["data"];
    let samples = this.datasets[this.activeDatasetIdx]["samples"];

    let samplesSettings = applyToolboxSettings(samples);

    // Rename and filter samples:
    this.filteredSettings = samplesSettings.filter((s) => !s.hidden);
    samples = this.filteredSettings.map((s) => s.name);
    data = data.filter((_, si) => !samplesSettings[si].hidden);

    return [data, samples];
  }

  resize(newHeight) {
    this.layout.height = newHeight;

    const maxTicks = (this.layout.height - 140) / 12;
    this.recalculateTicks(this.filteredSettings, this.layout.yaxis, maxTicks);

    super.resize(newHeight);
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let [data, samples] = this.prepData();
    if (data.length === 0 || samples.length === 0) return [];

    const maxTicks = (this.layout.height - 140) / 12;
    this.recalculateTicks(this.filteredSettings, this.layout.yaxis, maxTicks);

    let highlighted = this.filteredSettings.filter((s) => s.highlight);
    let traceParams = this.datasets[this.activeDatasetIdx]["trace_params"];

    return this.filteredSettings.map((sample, sampleIdx) => {
      let params = JSON.parse(JSON.stringify(traceParams)); // deep copy

      if (highlighted.length > 0) {
        if (sample.highlight !== null) {
          params.marker.color = sample.highlight;
        } else {
          params.marker.color = "grey";
        }
      }

      let values = data[sampleIdx];
      return {
        type: "box",
        x: values,
        name: sample.name,
        ...params,
      };
    });
  }

  exportData(format) {
    let [data, samples] = this.prepData();

    let delim = format === "tsv" ? "\t" : ",";

    let csv = "";
    for (let i = 0; i < data.length; i++) {
      csv += samples[i] + delim + data[i].join(delim) + "\n";
    }
    return csv;
  }
}
