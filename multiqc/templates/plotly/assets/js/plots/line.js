class LinePlot extends Plot {
  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx].lines.length; // no samples in a dataset
  }

  prepData() {
    // Prepare data to either build Plotly traces or export as a file
    let dataset = this.datasets[this.activeDatasetIdx];

    let lines = dataset.lines;

    let samples = lines.map((line) => line.name);
    let sampleSettings = applyToolboxSettings(samples);

    lines = lines.filter((line, idx) => {
      line.name = sampleSettings[idx].name ?? line.name;
      line.highlight = sampleSettings[idx].highlight;
      return !sampleSettings[idx].hidden;
    });

    return [samples, lines];
  }

  buildTraces(layout) {
    let dataset = this.datasets[this.activeDatasetIdx];

    let [samples, lines] = this.prepData();
    if (lines.length === 0 || samples.length === 0) return [];

    // Reorder points so highlighted points are on top
    let highlighted = lines.filter((p) => p.highlight);
    let nonHighlighted = lines.filter((p) => !p.highlight);
    lines = nonHighlighted.concat(highlighted);

    return lines.map((line) => {
      let x, y;
      if (line.data.length > 0 && Array.isArray(line.data[0])) {
        x = line.data.map((x) => x[0]);
        y = line.data.map((x) => x[1]);
      } else if (dataset.categories) {
        x = dataset.categories;
        y = line.data;
      } else {
        x = [...Array(line.data.length).keys()];
        y = line.data;
      }

      let params = JSON.parse(JSON.stringify(dataset["trace_params"])); // deep copy
      if (highlighted.length > 0) params.marker.color = line.highlight ?? "#cccccc";
      else params.marker.color = line.color;

      if (line["dashStyle"] !== undefined) params.line.dash = line["dashStyle"].toLowerCase();
      if (line["lineWidth"] !== undefined) params.line.width = line["lineWidth"];

      return {
        type: "scatter",
        x: x,
        y: y,
        name: line.name,
        text: x.map(() => line.name),
        ...params,
      };
    });
  }

  exportData(format) {
    let dataset = this.datasets[this.activeDatasetIdx];

    let [samples, lines] = this.prepData();

    // check if all lines have the same x values
    let sharedX = true;
    let x = null;
    lines.forEach((line) => {
      let thisX;
      if (line.data.length > 0 && Array.isArray(line.data[0])) {
        thisX = line.data.map((x) => x[0]);
      } else {
        thisX = dataset.categories;
      }
      if (x === null) {
        x = thisX;
      } else if (x.length !== thisX.length) {
        sharedX = false;
      } else if (x.some((v, i) => v !== thisX[i])) {
        sharedX = false;
      }
    });

    let sep = format === "tsv" ? "\t" : ",";
    let csv = "";
    if (sharedX) {
      csv += "Sample" + sep + x.join(sep) + "\n";
      lines.forEach((line) => {
        csv += line.name + sep + line.data.map((x) => x[1]).join(sep) + "\n";
      });
    } else {
      lines.forEach((line) => {
        csv += line.name + sep + "X" + sep + line.data.map((x) => x[0]).join(sep) + "\n";
        csv += line.name + sep + "Y" + sep + line.data.map((x) => x[1]).join(sep) + "\n";
      });
    }
    return csv;
  }
}