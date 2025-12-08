class LinePlot extends Plot {
  constructor(dump) {
    super(dump);
    this.heatmapMode = false;
    this.heatmapSortMode = "original"; // "original", "alphabetical", "clustered"
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx].lines.length; // no samples in a dataset
  }

  prepData(dataset) {
    // Prepare data to either build Plotly traces or export as a file
    dataset = dataset ?? this.datasets[this.activeDatasetIdx];

    let lines = dataset.lines;

    let samples = lines.map((line) => line.name);
    let sampleSettings = applyToolboxSettings(samples);
    this.filtSampleSettings = sampleSettings.filter((s) => !s.hidden);

    lines = lines.filter((line, idx) => {
      return !sampleSettings[idx].hidden;
    });

    lines = lines.map((line, idx) => {
      line.highlight = sampleSettings[idx].highlight;
      line.pseudonym = sampleSettings[idx].pseudonym;
      return line;
    });

    return [samples, lines];
  }

  plotAiHeader() {
    let result = super.plotAiHeader();
    if (this.pconfig.xlab) result += `X axis: ${this.pconfig.xlab}\n`;
    if (this.pconfig.ylab) result += `Y axis: ${this.pconfig.ylab}\n`;
    return result;
  }

  formatDatasetForAiPrompt(dataset) {
    let [samples, lines] = this.prepData(dataset);

    // Check if all samples are hidden
    if (samples.length === 0) {
      return "All samples are hidden by user, so no data to analyse. Please inform user to use the toolbox to unhide samples.\n";
    }

    const xsuffix = this.layout.xaxis.ticksuffix || "";
    const ysuffix = this.layout.yaxis.ticksuffix || "";

    let result = "Samples: " + samples.join(", ") + "\n\n";

    // If all y-values have the same suffix (like %), mention it in the header
    if (ysuffix) {
      result += `Y values are in ${ysuffix}\n\n`;
    }

    if (xsuffix) {
      result += `X values are in ${xsuffix}\n\n`;
    }

    const formattedLines = lines.map((line) => {
      let name = line.pseudonym ?? line.name;
      return {
        name: name,
        pairs: line.pairs.map((p) =>
          p.map((x, i) => {
            return !Number.isFinite(x) ? "" : Number.isInteger(x) ? x : parseFloat(x.toFixed(2));
          }),
        ),
      };
    });

    return (
      result +
      "\n\n" +
      formattedLines
        .map((line) => {
          return line.name + " " + line.pairs.map((p) => p.join(": ")).join(", ");
        })
        .join("\n\n")
    );
  }

  buildTraces() {
    if (this.heatmapMode) {
      return this.buildHeatmapTraces();
    }
    return this.buildLineTraces();
  }

  buildLineTraces() {
    let dataset = this.datasets[this.activeDatasetIdx];

    let [samples, lines] = this.prepData();
    if (lines.length === 0 || samples.length === 0) return [];

    // Reorder points so highlighted points are on top
    let highlighted = lines.filter((p) => p.highlight);
    let nonHighlighted = lines.filter((p) => !p.highlight);
    lines = nonHighlighted.concat(highlighted);

    return lines.map((line) => {
      let color = line.color;
      if (highlighted.length > 0) {
        color = line.highlight ?? "#cccccc";
      }

      let params = {
        line: {
          color: color,
          dash: line.dash,
          width: line.width,
        },
        marker: {
          color: color,
        },
        showlegend: line["showlegend"] ?? null,
        mode: line["mode"] ?? null,
      };

      let marker = line["marker"] ?? null;
      if (marker) {
        params.mode = "lines+markers";
        params.marker = {
          symbol: marker["symbol"],
          color: marker["fill_color"] ?? marker["color"] ?? color,
          line: {
            width: marker["width"],
            color: marker["line_color"] ?? marker["color"] ?? color,
          },
        };
      }

      updateObject(params, dataset["trace_params"], true);

      return {
        type: "scatter",
        x: line.pairs.map((x) => x[0]),
        y: line.pairs.map((x) => x[1]),
        name: line.name,
        text: line.pairs.map(() => line.name),
        ...params,
      };
    });
  }

  // Convert line data to heatmap format
  prepHeatmapData() {
    let [samples, lines] = this.prepData();
    if (lines.length === 0) return { rows: [], xcats: [], ycats: [] };

    // Get all unique x values across all lines, sorted
    let xSet = new Set();
    lines.forEach((line) => {
      line.pairs.forEach((pair) => xSet.add(pair[0]));
    });
    let xcats = Array.from(xSet);
    // Sort x categories (handle both numeric and string)
    xcats.sort((a, b) => {
      if (typeof a === "number" && typeof b === "number") return a - b;
      return String(a).localeCompare(String(b));
    });

    // Build a map from x to column index
    let xToIdx = new Map();
    xcats.forEach((x, i) => (xToIdx[x] = i));

    // Build rows: each line becomes a row
    let ycats = lines.map((line) => line.name);
    let rows = lines.map((line) => {
      let row = new Array(xcats.length).fill(null);
      line.pairs.forEach((pair) => {
        let colIdx = xToIdx[pair[0]];
        if (colIdx !== undefined) {
          row[colIdx] = pair[1];
        }
      });
      return row;
    });

    // Apply sorting
    if (this.heatmapSortMode === "alphabetical") {
      let indices = ycats.map((_, i) => i);
      indices.sort((a, b) => ycats[a].localeCompare(ycats[b]));
      ycats = indices.map((i) => ycats[i]);
      rows = indices.map((i) => rows[i]);
    } else if (this.heatmapSortMode === "clustered") {
      let order = this.clusterRows(rows);
      ycats = order.map((i) => ycats[i]);
      rows = order.map((i) => rows[i]);
    }

    return { rows, xcats, ycats };
  }

  // Simple hierarchical clustering (agglomerative, average linkage)
  clusterRows(rows) {
    if (rows.length <= 1) return rows.map((_, i) => i);

    // Compute distance matrix (Euclidean distance, handling nulls)
    const dist = (a, b) => {
      let sum = 0,
        count = 0;
      for (let i = 0; i < a.length; i++) {
        if (a[i] !== null && b[i] !== null) {
          sum += (a[i] - b[i]) ** 2;
          count++;
        }
      }
      return count > 0 ? Math.sqrt(sum / count) : Infinity;
    };

    // Initialize: each row is its own cluster
    let clusters = rows.map((row, i) => ({ indices: [i], centroid: row }));

    // Agglomerative clustering until we have one cluster
    while (clusters.length > 1) {
      // Find closest pair
      let minDist = Infinity,
        minI = 0,
        minJ = 1;
      for (let i = 0; i < clusters.length; i++) {
        for (let j = i + 1; j < clusters.length; j++) {
          let d = dist(clusters[i].centroid, clusters[j].centroid);
          if (d < minDist) {
            minDist = d;
            minI = i;
            minJ = j;
          }
        }
      }

      // Merge clusters minI and minJ
      let merged = {
        indices: clusters[minI].indices.concat(clusters[minJ].indices),
        centroid: clusters[minI].centroid.map((v, k) => {
          let a = clusters[minI].centroid[k];
          let b = clusters[minJ].centroid[k];
          if (a === null) return b;
          if (b === null) return a;
          return (a + b) / 2;
        }),
      };
      clusters.splice(minJ, 1);
      clusters.splice(minI, 1, merged);
    }

    return clusters[0].indices;
  }

  buildHeatmapTraces() {
    let { rows, xcats, ycats } = this.prepHeatmapData();
    if (rows.length === 0) return [];

    // Use Blues colorscale (monotonically increasing)
    let colorscale = [
      [0, "#f7fbff"],
      [0.125, "#deebf7"],
      [0.25, "#c6dbef"],
      [0.375, "#9ecae1"],
      [0.5, "#6baed6"],
      [0.625, "#4292c6"],
      [0.75, "#2171b5"],
      [0.875, "#08519c"],
      [1, "#08306b"],
    ];

    return [
      {
        type: "heatmap",
        z: rows,
        x: xcats,
        y: ycats,
        colorscale: colorscale,
        showscale: true,
        hoverongaps: false,
        hovertemplate: "<b>%{y}</b><br>%{x}: %{z}<extra></extra>",
      },
    ];
  }

  // Override afterPlotCreated to adjust layout for heatmap mode
  afterPlotCreated() {
    if (this.heatmapMode) {
      // Adjust y-axis for heatmap (categorical)
      Plotly.relayout(this.anchor, {
        "yaxis.type": "category",
        "yaxis.autorange": "reversed", // Put first sample at top
      });
    }
  }

  exportData(format) {
    let [_, lines] = this.prepData();

    // check if all lines have the same x values
    let sharedX = true;
    let x = null;
    lines.forEach((line) => {
      let thisX;
      thisX = line.pairs.map((x) => x[0]);
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
        csv += line.name + sep + line.pairs.map((x) => x[1]).join(sep) + "\n";
      });
    } else {
      lines.forEach((line) => {
        csv += line.name + sep + "X" + sep + line.pairs.map((x) => x[0]).join(sep) + "\n";
        csv += line.name + sep + "Y" + sep + line.pairs.map((x) => x[1]).join(sep) + "\n";
      });
    }
    return csv;
  }
}

// Make LinePlot globally available
window.LinePlot = LinePlot;
