class BoxPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.filteredSettings = [];
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx]["samples"].length;
  }

  prepData(dataset) {
    dataset = dataset ?? this.datasets[this.activeDatasetIdx];
    let data = dataset["data"];
    let samples = dataset["samples"];

    let samplesSettings = applyToolboxSettings(samples);

    // Rename and filter samples:
    this.filteredSettings = samplesSettings.filter((s) => !s.hidden);
    samples = this.filteredSettings.map((s) => s.name);
    data = data.filter((_, si) => !samplesSettings[si].hidden);

    return [data, samples];
  }

  formatDatasetForAiPrompt(dataset) {
    let prompt = "";

    // Prepare data to be sent to the LLM. LLM doesn't need things like colors, etc.
    let [data, samples] = this.prepData(dataset);

    // Check if all samples are hidden
    if (samples.length === 0) {
      return "All samples are hidden by user, so no data to analyse. Please inform user to use the toolbox to unhide samples.\n";
    }

    prompt += "|Sample|Min|Q1|Median|Q3|Max|Mean|\n";
    prompt += "|---|---|---|---|---|---|---|\n";

    const suffix = this.layout.xaxis.ticksuffix ? " " + this.layout.xaxis.ticksuffix : "";

    // Calculate statistics for each sample
    samples.forEach((sample, idx) => {
      const values = data[idx].filter((v) => Number.isFinite(v)).sort((a, b) => a - b);
      if (values.length === 0) return;

      let n = values.length;

      let min = values[0];
      let max = values[n - 1];
      let median = n % 2 === 1 ? values[Math.floor(n / 2)] : (values[n / 2 - 1] + values[n / 2]) / 2;
      let q1 = n >= 4 ? values[Math.floor(n / 4)] : values[0];
      let q3 = n >= 4 ? values[Math.floor((3 * n) / 4)] : values[n - 1];
      let mean = values.reduce((a, b) => a + b, 0) / n;

      // Format each value and add suffix
      let fmt = (val) => {
        if (!Number.isFinite(val)) return "";
        const isInt = Number.isInteger(val);
        return (isInt ? val : parseFloat(val.toFixed(2))) + suffix;
      };

      prompt += `|${anonymizeSampleName(sample)}|${fmt(min)}|${fmt(q1)}|${fmt(median)}|${fmt(q3)}|${fmt(max)}|${fmt(
        mean,
      )}|\n`;
    });

    return prompt;
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
