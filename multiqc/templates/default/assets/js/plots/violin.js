class ViolinPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.violinHeight = dump["violin_height"];
    this.extraHeight = dump["extra_height"];
    this.showTableByDefault = dump["show_table_by_default"];
    this.tableAnchor = dump["table_anchor"];
    this.isDownsampled = dump["is_downsampled"];
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx]["all_samples"].length;
  }

  prepData(dataset, keepHidden = false) {
    dataset = dataset ?? this.datasets[this.activeDatasetIdx];
    let metrics = dataset["metrics"];
    let headerByMetric = dataset["header_by_metric"];

    // Hidden metrics - check both the header and table configuration
    metrics = metrics.filter((metric) => {
      let header = headerByMetric[metric];
      // Check if column is hidden in table configuration
      let tableCheckbox = $(`#${this.tableAnchor}_config_modal_table .mqc_table_col_visible[value="${metric}"]`);
      let hiddenInTable = tableCheckbox.length > 0 && !tableCheckbox.is(":checked");

      return (header["hidden"] !== true || keepHidden) && !hiddenInTable;
    });

    let violinValuesBySampleByMetric = dataset["violin_value_by_sample_by_metric"];
    let scatterValuesBySampleByMetric = {};
    metrics.forEach((metric) => {
      let header = headerByMetric[metric];
      let scatterValuesBySample = {};
      if (header["show_points"]) {
        if (header["show_only_outliers"]) scatterValuesBySample = dataset["scatter_value_by_sample_by_metric"][metric];
        else scatterValuesBySample = violinValuesBySampleByMetric[metric];
      }
      scatterValuesBySampleByMetric[metric] = scatterValuesBySample;
    });

    let allSamples = dataset["all_samples"];
    let sampleSettings = applyToolboxSettings(allSamples);

    // Hidden samples
    let someHidden = sampleSettings.filter((s) => s.hidden).length > 0;
    if (someHidden) {
      let filteredViolinValuesBySampleByMetric = {};
      let filteredScatterValuesBySampleByMetric = {};
      metrics.map((metric) => {
        filteredViolinValuesBySampleByMetric[metric] = {};
        Object.keys(violinValuesBySampleByMetric[metric]).map((sample) => {
          if (!sampleSettings[allSamples.indexOf(sample)].hidden)
            filteredViolinValuesBySampleByMetric[metric][sample] = violinValuesBySampleByMetric[metric][sample];
        });
      });
      violinValuesBySampleByMetric = filteredViolinValuesBySampleByMetric;
      metrics.forEach((metric) => {
        let header = headerByMetric[metric];
        let scatterValuesBySample;
        if (header["show_points"] && header["show_only_outliers"]) {
          scatterValuesBySample = filteredScatterValuesBySampleByMetric[metric];
        } else {
          scatterValuesBySample = violinValuesBySampleByMetric[metric];
        }
        filteredScatterValuesBySampleByMetric[metric] = {};
        Object.keys(scatterValuesBySample).map((sample) => {
          if (!sampleSettings[allSamples.indexOf(sample)].hidden)
            filteredScatterValuesBySampleByMetric[metric][sample] = scatterValuesBySample[sample];
        });
      });
    }

    return [
      metrics,
      headerByMetric,
      allSamples,
      sampleSettings,
      violinValuesBySampleByMetric,
      scatterValuesBySampleByMetric,
    ];
  }

  plotAiHeader(view) {
    let prompt = "";
    if (view === "table") prompt += "Plot type: table\n";
    else prompt += "Plot type: violin plot\n";
    return prompt;
  }

  formatDatasetForAiPrompt(dataset) {
    let [
      metrics,
      headerByMetric,
      allSamples,
      sampleSettings,
      violinValuesBySampleByMetric,
      scatterValuesBySampleByMetric,
    ] = this.prepData(dataset, true);

    let results = "Number of samples: " + allSamples.length + "\n";
    if (this.isDownsampled)
      results +=
        "Note: sample number " +
        allSamples.length +
        " is greater than the threshold  so data points were downsampled to fit the context window. However, outliers for each metric were identified and kept in the datasets.\n";
    results += "\n";

    if (metrics.length === 0) {
      results +=
        'All columns are hidden by user, so no data to analyse. Please inform user to use the "Configure columns" button to make some columns visible.\n';
      return results;
    }

    // Check if all samples are hidden
    if (sampleSettings.every((s) => s.hidden)) {
      results +=
        "All samples are hidden by user, so no data to analyse. Please inform user to use the toolbox to unhide samples.\n";
      return results;
    }

    results += "Metrics:\n";
    results += metrics
      .map((metric) => `${headerByMetric[metric].title} - ${headerByMetric[metric].description}`)
      .join("\n");
    results += "\n\n";

    results +=
      `|${this.pconfig.col1_header}|` + metrics.map((metric) => headerByMetric[metric].title).join("|") + "|\n";
    results += "|---|" + metrics.map(() => "---").join("|") + "|\n";
    results += sampleSettings
      .map((sample) => {
        if (sample.hidden) return "";
        if (
          metrics.every(
            (metric) =>
              violinValuesBySampleByMetric[metric][sample.originalName] === undefined &&
              scatterValuesBySampleByMetric[metric][sample.originalName] === undefined,
          )
        )
          return "";

        return (
          `|${sample.pseudonym ?? sample.name}|` +
          metrics
            .map((metric) => {
              const value =
                violinValuesBySampleByMetric[metric][sample.originalName] ??
                scatterValuesBySampleByMetric[metric][sample.originalName];
              const suffix = headerByMetric[metric].suffix;
              if (value === undefined || value === null) return "";
              if (typeof value === "string") return value + (suffix ?? "");
              if (Number.isFinite(value)) {
                if (Number.isInteger(value)) return value;
                return value.toFixed(2);
              }
              return "";
            })
            .join("|") +
          "|\n"
        );
      })
      .join("");
    return results;
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let dataset = this.datasets[this.activeDatasetIdx];
    let [
      metrics,
      headerByMetric,
      allSamples,
      sampleSettings,
      violinValuesBySampleByMetric,
      scatterValuesBySampleByMetric,
    ] = this.prepData();

    let outliersWarning = true;
    metrics.forEach((metric) => {
      let header = headerByMetric[metric];
      if (!header["show_points"] || !header["show_only_outliers"]) {
        outliersWarning = false;
      }
    });
    if (outliersWarning)
      $("#table-violin-info-" + this.anchor).append(" For efficiency, separate points are shown only for outliers.");

    let layout = this.layout;
    layout.height = this.violinHeight * metrics.length + this.extraHeight;
    $("#" + this.anchor + "-wrapper").css("height", layout.height + "px");
    if (metrics.length === 0) return [];

    layout.grid.rows = metrics.length;
    layout.grid.subplots = metrics.map((metric, metricIdx) => {
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;
      return ["x" + axisKey + "y" + axisKey];
    });

    metrics.map((metric, metricIdx) => {
      let header = headerByMetric[metric];

      // Set layouts for each violin individually
      layout["yaxis" + (metricIdx + 1)] = {
        automargin: layout["yaxis"]["automargin"],
        color: layout["yaxis"]["color"],
        gridcolor: layout["yaxis"]["gridcolor"],
        zerolinecolor: layout["yaxis"]["zerolinecolor"],
        hoverformat: layout["yaxis"]["hoverformat"],
        tickfont: {
          size: layout["yaxis"]["tickfont"]["size"],
          color: layout["yaxis"]["tickfont"]["color"],
        },
      };
      layout["xaxis" + (metricIdx + 1)] = {
        automargin: layout["xaxis"]["automargin"],
        color: layout["xaxis"]["color"],
        gridcolor: layout["xaxis"]["gridcolor"],
        zerolinecolor: layout["xaxis"]["zerolinecolor"],
        hoverformat: layout["xaxis"]["hoverformat"],
        tickfont: {
          size: layout["xaxis"]["tickfont"]["size"],
          color: layout["xaxis"]["tickfont"]["color"],
        },
      };
      if (header["xaxis"] !== undefined && header["xaxis"] !== null) {
        layout["xaxis" + (metricIdx + 1)] = Object.assign(layout["xaxis" + (metricIdx + 1)], header["xaxis"]);
      }
      let title = header.title + "  ";
      if (header["namespace"]) title = header["namespace"] + "  <br>" + title;
      layout["yaxis" + (metricIdx + 1)]["tickmode"] = "array";
      layout["yaxis" + (metricIdx + 1)]["tickvals"] = [metricIdx];
      layout["yaxis" + (metricIdx + 1)]["ticktext"] = [title];

      if (header["hoverformat"] !== undefined && header["hoverformat"] !== null) {
        layout["xaxis" + (metricIdx + 1)]["hoverformat"] = header["hoverformat"];
      }

      // Set color for each violin individually
      if (header["color"] !== undefined && header["color"] !== null) {
        layout["yaxis" + (metricIdx + 1)]["tickfont"] = { color: "rgb(" + header["color"] + ")" };
      }
    });

    layout.xaxis = layout["xaxis1"];
    layout.yaxis = layout["yaxis1"];

    let traces = [];
    metrics.map((metric, metricIdx) => {
      let header = headerByMetric[metric];
      let params = JSON.parse(JSON.stringify(dataset["trace_params"])); // deep copy

      // Set color for each violin individually
      if (header["color"]) {
        params["fillcolor"] = "rgb(" + header["color"] + ")";
        params["line"]["color"] = "rgb(" + header["color"] + ")";
      }

      // Create violin traces
      let violinValuesBySample = violinValuesBySampleByMetric[metric];
      let samples = [],
        values = [];
      Object.entries(violinValuesBySample).map(([sample, value]) => {
        samples.push(sample);
        values.push(value);
      });

      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;
      traces.push({
        type: "violin",
        x: values,
        name: metricIdx, // headerByMetric[metric].title + "  ",
        text: samples, // sample names
        xaxis: "x" + axisKey,
        yaxis: "y" + axisKey,
        ...params,
      });
    });

    // We want to select points on all violins belonging to this specific sample.
    // Points are rendered each as a separate trace, so we need to collect
    // each trace (curve) number for each point by sample.
    let currentCurveNumber = traces.length;
    let curveNumbersBySample = {};
    let curveAxisBySample = {};
    let scatterDataByMetric = [];
    metrics.map((metric, metricIdx) => {
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;
      let valuesBySample = scatterValuesBySampleByMetric[metric];

      let scatterData = [];
      Object.entries(valuesBySample).map(([sample, value]) => {
        scatterData.push([sample, value]);

        if (!curveNumbersBySample[sample]) {
          curveNumbersBySample[sample] = [];
          curveAxisBySample[sample] = [];
        }
        curveNumbersBySample[sample].push(currentCurveNumber++);
        curveAxisBySample[sample].push("x" + axisKey + "y" + axisKey);
      });
      scatterDataByMetric.push(scatterData);
    });

    let highlightingEnabled = sampleSettings.filter((s) => s.highlight).length > 0;

    let seed = 42.1231;
    function random() {
      // Math.random does not have a seed, so we use this
      let x = Math.sin(seed++) * 10000;
      return x - Math.round(x);
    }
    let scatters = [];
    // Add scatter plots on top of violins to show individual points
    // Plotly supports showing points automatically with `points="all"`,
    // however, it's problematic to give each sample individual color,
    // and set up hover events to highlight all points for a sample. So as
    // a workaround, we add a separate scatter plot. One thing to solve later
    // would be to be able to only show outliers, not all points, as the violin
    // plot can do that automatically, and we lose this functionality here.
    scatterDataByMetric.map((scatterData, metricIdx) => {
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;

      scatterData.map(([sample, value]) => {
        let state = sampleSettings[allSamples.indexOf(sample)];
        let params = JSON.parse(JSON.stringify(dataset["scatter_trace_params"])); // deep copy

        let color = params.marker.color;
        let size = params.marker.size;
        if (highlightingEnabled) {
          color = state.highlight ?? "#ddd";
          size = state.highlight !== null ? 8 : size;
        }

        let customData = {
          curveNumbers: curveNumbersBySample[sample],
          curveAxis: curveAxisBySample[sample],
        };

        const jitter = 0.3;
        scatters.push({
          type: "scatter",
          x: [value],
          y: [metricIdx + random() * jitter], // add vertical jitter
          text: [state.name ?? sample],
          xaxis: "x" + axisKey,
          yaxis: "y" + axisKey,
          customdata: customData,
          ...params,
          marker: {
            color: color,
            size: size,
          },
        });
      });
    });
    traces = traces.concat(scatters);
    return traces;
  }

  exportData(format) {
    let [
      metrics,
      headerByMetric,
      allSamples,
      sampleSettings,
      violinValuesBySampleByMetric,
      scatterValuesBySampleByMetric,
    ] = this.prepData();

    let sep = format === "tsv" ? "\t" : ",";
    // Export all data points as a table, samples are rows, metrics are columns
    let titles = metrics.map((metric) => headerByMetric[metric].title);
    // Escape titles that contain the separator character
    titles = titles.map((title) => (title.includes(sep) ? `"${title}"` : title));
    let csv = "Sample" + sep + titles.join(sep) + "\n";
    for (let i = 0; i < allSamples.length; i++) {
      let sample = allSamples[i];
      if (sampleSettings[i].hidden) continue;
      csv += sample + sep;
      csv += metrics
        .map((metric) => {
          let val = violinValuesBySampleByMetric[metric][sample];
          if (val === undefined) val = ".";
          return val;
        })
        .join(sep);
      csv += "\n";
    }
    return csv;
  }

  afterPlotCreated() {
    let anchor = this.anchor;
    let plot = document.getElementById(anchor);

    plot.on("plotly_hover", function (eventdata) {
      if (!eventdata.points) return;
      let point = eventdata.points[0];
      if (point.data.type === "scatter") {
        let curveNumbers = point.data.customdata["curveNumbers"];
        let curveAxis = point.data.customdata["curveAxis"];
        let points = curveNumbers.map((curveNum) => {
          return { curveNumber: curveNum, pointNumber: 0 };
        });
        // console.log("hover", point.curveNumber);
        // hoverInfoDiv.html(point.data.text);
        Plotly.Fx.hover(anchor, points, curveAxis);
      }
    });
  }
}
