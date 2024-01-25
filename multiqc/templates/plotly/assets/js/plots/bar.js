class BarPlot extends Plot {
  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    let cats = this.datasets[this.activeDatasetIdx]["cats"];
    if (cats.length === 0) return 0; // no categories
    return cats[0].data.length; // no data for a category
  }

  prepData() {
    let cats = this.datasets[this.activeDatasetIdx]["cats"];
    let samples = this.datasets[this.activeDatasetIdx]["samples"];

    let samplesSettings = applyToolboxSettings(samples);

    // Rename and filter samples:
    let filteredSettings = samplesSettings.filter((s) => !s.hidden);
    samples = filteredSettings.map((s) => s.name);

    cats = cats.map((cat) => {
      let data = this.pActive ? cat["data_pct"] : cat.data;
      return {
        data: data.filter((_, si) => !samplesSettings[si].hidden),
        color: cat.color,
        name: cat.name,
      };
    });

    return [samples, cats, filteredSettings];
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let [samples, cats, filteredSettings] = this.prepData();
    if (cats.length === 0 || samples.length === 0) return [];

    return cats.map((cat) => {
      let params = JSON.parse(JSON.stringify(this.traceParams)); // deep copy
      params.marker.color = cat.color;
      params.marker.line = {
        // Remove grey from highlights:
        color: filteredSettings.map((s) => (s.highlight && s.highlight !== "#cccccc" ? s.highlight : null)),
        width: filteredSettings.map((s) => (s.highlight && s.highlight !== "#cccccc" ? 2 : params.marker.line.width)),
      };

      return {
        type: "bar",
        x: cat.data,
        y: samples,
        name: cat.name,
        ...params,
      };
    });
  }

  exportData(format) {
    let [samples, cats, _] = this.prepData();

    let delim = format === "tsv" ? "\t" : ",";

    let csv = "Sample" + delim + cats.map((cat) => cat.name).join(delim) + "\n";
    for (let i = 0; i < samples.length; i++) {
      csv += samples[i] + delim + cats.map((cat) => cat.data[i]).join(delim) + "\n";
    }
    return csv;
  }
}
