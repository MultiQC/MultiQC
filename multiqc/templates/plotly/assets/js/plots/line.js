class LinePlot extends Plot {
  constructor(dump) {
    super(dump);

    this.categories = dump["categories"];
    this.yAutorangeBeforeBands = dump["y_autorange_before_bands"];

    // Tracking Y-axis range to maintain the "Y-Limits" toggle button
    this.ymin = this.layout.yaxis.range[0];
    this.ymax = this.layout.yaxis.range[1];
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx].lines.length; // no samples in a dataset
  }

  prepData() {
    let lines = this.datasets[this.activeDatasetIdx].lines;

    let samples = lines.map((line) => line.name);
    let sampleSettings = applyToolboxSettings(samples);
    if (sampleSettings == null) return []; // All series are hidden, do not render the graph.
    lines = lines.filter((line, idx) => {
      line.name = sampleSettings[idx].name ?? line.name;
      line.highlight = sampleSettings[idx].highlight;
      return !sampleSettings[idx].hidden;
    });

    return [samples, lines];
  }

  buildTraces() {
    let [samples, lines] = this.prepData();
    if (lines.length === 0 || samples.length === 0) return [];

    // Reorder points so highlighted points are on top
    let highlighted = lines.filter((p) => p.highlight);
    let nonHighlighted = lines.filter((p) => !p.highlight);
    lines = nonHighlighted.concat(highlighted);

    // Toggle buttons for Y-axis limis
    let ymaxSet = this.ymax !== "undefined" && this.ymax !== null;
    let yminSet = this.ymin !== "undefined" && this.ymin !== null;
    // Only create if there is a y-axis limit
    let groupDiv = $("#" + this.target).closest(".mqc_hcplot_plotgroup");
    let addLimitToggle = ymaxSet || (yminSet && this.ymin !== 0);
    if (this.yAutorangeBeforeBands) {
      if (
        yminSet &&
        this.yAutorangeBeforeBands[0] <= this.ymin &&
        ymaxSet &&
        this.yAutorangeBeforeBands[1] >= this.ymax
      )
        addLimitToggle = false;
    }
    if (addLimitToggle) {
      let wrapper = $('<div class="mqc_hcplot_yaxis_limit_toggle hidden-xs" />').prependTo(groupDiv);
      wrapper.append(
        '<span class="mqc_switch_wrapper"' +
          '"' +
          ' data-ymax="' +
          this.ymax +
          '"' +
          ' data-ymin="' +
          this.ymin +
          '"' +
          ' data-target="' +
          this.target +
          '"' +
          ' data-y_autorange_range_before_bands="' +
          this.yAutorangeBeforeBands +
          '">Y-Limits: <span class="mqc_switch on">on</span></span>',
      );
      wrapper.after('<div class="clearfix" />');
    }

    return lines.map((line) => {
      let x, y;
      if (line.data.length > 0 && Array.isArray(line.data[0])) {
        x = line.data.map((x) => x[0]);
        y = line.data.map((x) => x[1]);
      } else {
        x = [...Array(line.data.length).keys()];
        y = line.data;
      }

      let params = JSON.parse(JSON.stringify(this.traceParams)); // deep copy
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
    let [samples, lines] = this.prepData();

    // check if all lines have the same x values
    let sameX = true;
    let x = null;
    lines.forEach((line) => {
      let thisX;
      if (line.data.length > 0 && Array.isArray(line.data[0])) {
        thisX = line.data.map((x) => x[0]);
      } else {
        thisX = this.categories;
      }
      if (x === null) {
        x = thisX;
      } else if (x.length !== thisX.length) {
        sameX = false;
      } else if (x.some((v, i) => v !== thisX[i])) {
        sameX = false;
      }
    });

    let delim = format === "tsv" ? "\t" : ",";
    let csv = "";
    if (sameX) {
      // can create a CSV where columns are samples, and rows are x values, and cells are y values
      csv += "Sample" + delim + lines.map((line) => line.name).join(delim) + "\n";
      for (let i = 0; i < x.length; i++) {
        csv += x[i] + delim + lines.map((line) => line.data[i][1]).join(delim) + "\n";
      }
    } else {
      // can create a CSV where columns are: sample - x - y
      csv += "Sample" + delim + "X" + delim + "Y" + "\n";
      lines.forEach((line) => {
        csv +=
          line.name +
          delim +
          line.data.map((x) => x[0]).join(";") +
          delim +
          line.data.map((x) => x[1]).join(";") +
          "\n";
      });
    }
    return csv;
  }
}
