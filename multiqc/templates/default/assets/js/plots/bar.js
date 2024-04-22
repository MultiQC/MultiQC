class BarPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.filteredSettings = [];
  }

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
    this.filteredSettings = samplesSettings.filter((s) => !s.hidden);
    samples = this.filteredSettings.map((s) => s.name);

    cats = cats.map((cat) => {
      let data = this.pActive ? cat["data_pct"] : cat.data;
      return {
        data: data.filter((_, si) => !samplesSettings[si].hidden),
        color: cat.color, // formatted as "r,g,b", to be wrapped with "rgb()" or "rgba()"
        name: cat.name,
      };
    });

    return [cats, samples];
  }

  resize(newHeight) {
    this.layout.height = newHeight;

    const maxTicks = (this.layout.height - 140) / 12;
    this.recalculateTicks(this.filteredSettings, this.layout.yaxis, maxTicks);

    super.resize(newHeight);
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let [cats, samples] = this.prepData();
    if (cats.length === 0 || samples.length === 0) return [];

    const maxTicks = (this.layout.height - 140) / 12;
    this.recalculateTicks(this.filteredSettings, this.layout.yaxis, maxTicks);

    let highlighted = this.filteredSettings.filter((s) => s.highlight);
    let firstHighlightedSample = this.firstHighlightedSample(this.filteredSettings);
    let traceParams = this.datasets[this.activeDatasetIdx]["trace_params"];

    return cats.map((cat) => {
      if (this.layout.barmode !== "group") {
        // Plotting each sample as a separate trace to be able to set alpha for each
        // sample color separately, so we can dim the de-highlighted samples.
        return this.filteredSettings.map((sample, sampleIdx) => {
          let params = JSON.parse(JSON.stringify(traceParams)); // deep copy

          let alpha = highlighted.length > 0 && sample.highlight === null ? 0.1 : 1;
          params.marker.color = "rgba(" + cat.color + "," + alpha + ")";

          return {
            type: "bar",
            x: [cat.data[sampleIdx]],
            y: [sample.name],
            name: cat.name,
            meta: cat.name,
            // To make sure the legend uses bright category colors and not the dim ones:
            showlegend: sampleIdx === firstHighlightedSample,
            legendgroup: cat.name,
            ...params,
          };
        });
      } else {
        // "group"
        // Plotly adds giant gaps between bars in the group mode when adding each sample as a
        // separate trace. Sacrificing dimming the de-highlighted bars to get rid of this gap.
        let params = JSON.parse(JSON.stringify(traceParams)); // deep copy
        samples = this.filteredSettings.map((s) => s.name);
        params.marker.color = "rgb(" + cat.color + ")";

        return {
          type: "bar",
          x: cat.data,
          y: samples,
          name: cat.name,
          meta: cat.name,
          ...params,
        };
      }
    });
  }

  exportData(format) {
    let [cats, filteredSettings] = this.prepData();

    let delim = format === "tsv" ? "\t" : ",";

    let csv = "Sample" + delim + cats.map((cat) => cat.name).join(delim) + "\n";
    for (let i = 0; i < filteredSettings.length; i++) {
      csv += filteredSettings[i].name + delim + cats.map((cat) => cat.data[i]).join(delim) + "\n";
    }
    return csv;
  }
}
