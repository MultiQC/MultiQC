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
    let filteredSettings = samplesSettings.filter((s) => !s.hidden);

    cats = cats.map((cat) => {
      let data = this.pActive ? cat["data_pct"] : cat.data;
      return {
        data: data.filter((_, si) => !samplesSettings[si].hidden),
        color: cat.color,
        name: cat.name,
      };
    });

    return [cats, filteredSettings];
  }

  resize(newHeight) {
    if (newHeight === null || newHeight === undefined) console.error("BarPlot.resize: newHeight is " + newHeight);

    this.layout.height = newHeight;
    this.recalculateTicks();
    super.resize(newHeight);
  }

  recalculateTicks() {
    function subsample(values, num, start = 0) {
      // Take ~`num` samples from values evenly, always include `start`
      if (values.length <= num) return values;

      let indices = Array.from({ length: values.length }, (_, i) => i);
      let after = indices.slice(start);
      let before = indices.slice(0, start + 1); // including the pivot
      before.reverse(); // will stepping back
      let step = Math.floor(values.length / num);
      before = before.filter((_, i) => i % step === 0);
      after = after.filter((_, i) => i % step === 0);
      before.reverse();
      before = before.slice(0, before.length - 1); // remove the pivot
      indices = before.concat(after);
      return indices.map((i) => values[i]);
    }

    let highlighted = this.filteredSettings.filter((s) => s.highlight);
    let firstHighlightedSample = this.firstHighlightedSample();

    if (highlighted.length === 0) {
      this.layout.yaxis.tickmode = null;
      this.layout.yaxis.tickvals = null;
      this.layout.yaxis.ticktext = null;
    } else {
      // Have to switch to tickmode=array to set colors to ticks. however, this way plotly will try
      // to fit _all_ ticks on the screen, and if there are too many, they will overlap. to prevent that,
      // if there are too many samples, we will show only highlighted samples plus a subsampled number
      // of ticks, but up to a constant:
      this.layout.yaxis.tickmode = "array";
      if (this.layout.height === null || this.layout.height === undefined)
        console.error("BarPlot.recalculateTicks: this.layout.height is " + this.layout.height);

      const maxTicks = (this.layout.height - 140) / 15; // 20px per tick
      let selected = subsample(this.filteredSettings, maxTicks, firstHighlightedSample);

      this.layout.yaxis.tickvals = selected.map((s) => s.name);
      this.layout.yaxis.ticktext = selected.map((s) => "<span style='color:" + s.highlight + "'>" + s.name + "</span>");
    }
  }

  firstHighlightedSample() {
    let index = 0;
    let highlighted = this.filteredSettings.filter((s) => s.highlight);
    if (highlighted.length > 0) index = this.filteredSettings.findIndex((s) => s.highlight);
    return index;
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces(layout) {
    let [cats, filteredSettings] = this.prepData();
    if (cats.length === 0 || filteredSettings.length === 0) return [];
    this.filteredSettings = filteredSettings;

    this.recalculateTicks();

    let highlighted = filteredSettings.filter((s) => s.highlight);
    let firstHighlightedSample = this.firstHighlightedSample();
    let traceParams = this.datasets[this.activeDatasetIdx]["trace_params"];

    return cats.map((cat) => {
      return filteredSettings.map((sample, sampleIdx) => {
        let params = JSON.parse(JSON.stringify(traceParams)); // deep copy

        let alpha = highlighted.length > 0 && !sample.highlight ? 0.1 : 1;
        params.marker.color = "rgba(" + cat.color + "," + alpha + ")";

        return {
          type: "bar",
          x: [cat.data[sampleIdx]],
          y: [sample.name],
          name: cat.name,
          meta: cat.name,
          // To make sure the legend uses bright category colors and not the deemed ones.
          showlegend: sampleIdx === firstHighlightedSample,
          legendgroup: cat.name,
          ...params,
        };
      });
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
