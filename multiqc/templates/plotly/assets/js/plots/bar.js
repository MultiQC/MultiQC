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

  // Constructs and returns traces for the Plotly plot
  buildTraces(layout) {
    let [cats, filteredSettings] = this.prepData();
    if (cats.length === 0 || filteredSettings.length === 0) return [];

    // Use subplots only to set different color to each bar label when using
    // highlight interactivity
    layout.grid = {
      columns: 1,
      roworder: "top to bottom",
      ygap: 0,
      rows: filteredSettings.length,
      subplots: filteredSettings.map((_, sampleIdx) => {
        return ["xy" + (sampleIdx === 0 ? "" : sampleIdx + 1)];
      }),
    };

    let anyHighlight = filteredSettings.some((s) => s.highlight);

    filteredSettings.forEach((sample, sampleIdx) => {
      layout["yaxis" + (sampleIdx + 1)] = {
        showgrid: false,
        gridcolor: layout["yaxis"]["gridcolor"],
        zerolinecolor: layout["yaxis"]["zerolinecolor"],
        color: layout["yaxis"]["color"],
        tickfont: {
          size: layout["yaxis"]["tickfont"]["size"],
          color: anyHighlight ? sample.highlight ?? "#cccccc" : "black",
        },
        hoverformat: layout["yaxis"]["hoverformat"],
        ticksuffix: layout["yaxis"]["ticksuffix"],
        automargin: true,
      };
    });
    layout.yaxis = layout["yaxis1"];

    let traceParams = this.datasets[this.activeDatasetIdx]["trace_params"];

    // Because we are using subplots, we need to make sure the legend is not
    // duplicated for each sample, so we pick one sample for that as all
    // samples have data for all categories.
    let legendSample = 0;
    if (anyHighlight)
      // To make sure the legend uses the bright colors from the toolbox
      legendSample = filteredSettings.findIndex((s) => s.highlight);

    return cats.map((cat) => {
      return filteredSettings.map((sample, sampleIdx) => {
        let params = JSON.parse(JSON.stringify(traceParams)); // deep copy

        let alpha = anyHighlight && !sample.highlight ? 0.1 : 1;
        params.marker.color = "rgba(" + cat.color + "," + alpha + ")";

        return {
          type: "bar",
          x: [cat.data[sampleIdx]],
          y: [sample.name],
          name: cat.name,
          yaxis: "y" + (sampleIdx === 0 ? "" : sampleIdx + 1),
          ...params,
          legendgroup: cat.name,
          showlegend: sampleIdx === legendSample,
          meta: cat.name,
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