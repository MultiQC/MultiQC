class BoxPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.filteredSettings = [];
    this.sortSwitchSortedActive = dump["sort_switch_sorted_active"];
    this.isStatsData = dump.datasets && dump.datasets.length > 0 ? dump.datasets[0].is_stats_data : false;
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx]["samples"].length;
  }

  prepData(dataset) {
    dataset = dataset ?? this.datasets[this.activeDatasetIdx];

    // Choose data and samples based on sorting state
    let data = this.sortSwitchSortedActive && dataset["data_sorted"] ? dataset["data_sorted"] : dataset["data"];
    let samples =
      this.sortSwitchSortedActive && dataset["samples_sorted"] ? dataset["samples_sorted"] : dataset["samples"];

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

    // Format each value and add suffix
    let fmt = (val) => {
      if (!Number.isFinite(val)) return "";
      const isInt = Number.isInteger(val);
      return (isInt ? val : parseFloat(val.toFixed(2))) + suffix;
    };

    // Handle statistics or raw data
    samples.forEach((sample, idx) => {
      if (this.isStatsData) {
        // Use pre-calculated statistics
        const stats = data[idx];
        const min = stats.min || 0;
        const max = stats.max || 0;
        const median = stats.median || 0;
        const q1 = stats.q1 || min;
        const q3 = stats.q3 || max;
        const mean = stats.mean || median;

        prompt += `|${sample}|${fmt(min)}|${fmt(q1)}|${fmt(median)}|${fmt(q3)}|${fmt(max)}|${fmt(mean)}|\n`;
      } else {
        // Calculate statistics for raw data
        const values = data[idx].filter((v) => Number.isFinite(v)).sort((a, b) => a - b);
        if (values.length === 0) return;

        let n = values.length;

        let min = values[0];
        let max = values[n - 1];
        let median = n % 2 === 1 ? values[Math.floor(n / 2)] : (values[n / 2 - 1] + values[n / 2]) / 2;
        let q1 = n >= 4 ? values[Math.floor(n / 4)] : values[0];
        let q3 = n >= 4 ? values[Math.floor((3 * n) / 4)] : values[n - 1];
        let mean = values.reduce((a, b) => a + b, 0) / n;

        prompt += `|${sample}|${fmt(min)}|${fmt(q1)}|${fmt(median)}|${fmt(q3)}|${fmt(max)}|${fmt(mean)}|\n`;
      }
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

      if (this.isStatsData) {
        // Create box plot from statistics
        return {
          type: "box",
          q1: [values.q1 || 0],
          median: [values.median || 0],
          q3: [values.q3 || 0],
          lowerfence: [values.min || 0],
          upperfence: [values.max || 0],
          mean: [values.mean || values.median || 0],
          y: [sample.name], // Add y-coordinate for horizontal box plot positioning
          name: sample.name,
          ...params,
        };
      } else {
        // Create box plot from raw data
        return {
          type: "box",
          x: values,
          name: sample.name,
          ...params,
        };
      }
    });
  }

  exportData(format) {
    let [data, samples] = this.prepData();

    let delim = format === "tsv" ? "\t" : ",";

    let csv = "";

    if (this.isStatsData) {
      // Export statistics as CSV
      csv =
        "Sample" + delim + "Min" + delim + "Q1" + delim + "Median" + delim + "Q3" + delim + "Max" + delim + "Mean\n";
      for (let i = 0; i < data.length; i++) {
        const stats = data[i];
        csv +=
          samples[i] +
          delim +
          (stats.min || 0) +
          delim +
          (stats.q1 || 0) +
          delim +
          (stats.median || 0) +
          delim +
          (stats.q3 || 0) +
          delim +
          (stats.max || 0) +
          delim +
          (stats.mean || stats.median || 0) +
          "\n";
      }
    } else {
      // Export raw data
      for (let i = 0; i < data.length; i++) {
        csv += samples[i] + delim + data[i].join(delim) + "\n";
      }
    }

    return csv;
  }
}

// Make BoxPlot globally available
window.BoxPlot = BoxPlot;

$(function () {
  // Listener for box plot sorting toggle - use exact same pattern as heatmap
  $('button[data-action="unsorted"], button[data-action="sorted_by_median"]').on("click", function (e) {
    e.preventDefault();
    let $btn = $(this);
    let plotAnchor = $(this).data("plot-anchor");
    let plot = mqc_plots[plotAnchor];

    // Only proceed if this is a box plot
    if (!plot || !(plot instanceof BoxPlot)) {
      return;
    }

    // Toggle buttons
    $btn.toggleClass("active").siblings().toggleClass("active");

    // Update plot state
    plot.sortSwitchSortedActive = $btn.data("action") === "sorted_by_median";

    // Re-render plot
    renderPlot(plotAnchor);
  });
});
