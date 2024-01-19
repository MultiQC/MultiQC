class ViolinPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.scatterTraceParams = dump["scatter_trace_params"];
    this.showOnlyOutliers = dump["show_only_outliers"];
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx]["all_samples"].length;
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    if (this.showOnlyOutliers) {
      $("#table-violin-info-" + this.target).append(" For efficiency, separate points are shown only for outliers.");
    }

    let dataset = this.datasets[this.activeDatasetIdx];
    let metrics = dataset["metrics"];
    if (metrics.length === 0) return [];
    let headerByMetric = dataset["header_by_metric"];
    let valuesBySampleByMetric = dataset["values_by_sample_by_metric"];
    let outliersByMetric = dataset["outliers_by_metric"];
    let layout = this.layout;
    const traceParams = this.traceParams;
    const scatterTraceParams = this.scatterTraceParams;

    let allSamples = this.datasets[this.activeDatasetIdx]["all_samples"];
    let sampleSettings = applyToolboxSettings(allSamples);
    if (sampleSettings == null) return; // All series are hidden, do not render the graph.

    // Hidden metrics
    metrics = metrics.filter((metric) => {
      let header = headerByMetric[metric];
      return header["hidden"] !== true;
    });
    layout.height = 70 * metrics.length + 50;
    $("#" + this.target + "-wrapper").css("height", layout.height + "px");
    if (metrics.length === 0) {
      return [];
    }
    layout.grid.rows = metrics.length;
    layout.grid.subplots = metrics.map((metric, metricIdx) => {
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;
      return ["x" + axisKey + "y" + axisKey];
    });

    // Hidden samples
    let filteredValuesBySampleByMetric = {};
    metrics.map((metric) => {
      filteredValuesBySampleByMetric[metric] = {};
      Object.keys(valuesBySampleByMetric[metric]).map((sample, sampleIdx) => {
        if (!sampleSettings[allSamples.indexOf(sample)].hidden) {
          filteredValuesBySampleByMetric[metric][sample] = valuesBySampleByMetric[metric][sample];
        }
      });
    });
    valuesBySampleByMetric = filteredValuesBySampleByMetric;

    let traces = [];
    metrics.map((metric, metricIdx) => {
      let params = JSON.parse(JSON.stringify(traceParams)); // deep copy
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;

      let header = headerByMetric[metric];
      let valuesBySample = valuesBySampleByMetric[metric];
      let outliers = outliersByMetric[metric];

      // Set layouts for each violin individually
      let title = header.title + "  ";
      if (header["namespace"]) title = header["namespace"] + "  <br>" + title;
      layout["xaxis" + (metricIdx + 1)] = {
        gridcolor: layout["xaxis"]["gridcolor"],
        zerolinecolor: layout["xaxis"]["zerolinecolor"],
        tickfont: {
          size: 9,
          color: "rgba(0,0,0,0.5)",
        },
      };
      layout["yaxis" + (metricIdx + 1)] = {
        gridcolor: layout["yaxis"]["gridcolor"],
        zerolinecolor: layout["yaxis"]["zerolinecolor"],
        automargin: true,
      };
      if (header["xaxis"] !== undefined) {
        layout["xaxis" + (metricIdx + 1)] = Object.assign(layout["xaxis" + axisKey], header["xaxis"]);
      }
      layout["yaxis" + (metricIdx + 1)]["tickmode"] = "array";
      layout["yaxis" + (metricIdx + 1)]["tickvals"] = [metricIdx];
      layout["yaxis" + (metricIdx + 1)]["ticktext"] = [title];

      params["line"] = { width: 4 };

      // Set color for each violin individually
      if (header["color"] !== undefined) {
        layout["yaxis" + axisKey]["tickfont"] = { color: "rgba(" + header["color"] + ",1)" };
        params["fillcolor"] = "rgba(" + header["color"] + ",1)";
        params["line"]["color"] = "rgba(" + header["color"] + ",1)";
      }

      // Create violin traces
      let samples = [],
        values = [];
      Object.entries(valuesBySample).map(([sample, value]) => {
        samples.push(sample);
        values.push(value);
      });

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

    layout["xaxis"] = layout["xaxis1"];
    layout["yaxis"] = layout["yaxis1"];

    // We want to select points on all violins belonging to this specific sample.
    // Points are rendered each as a separate trace, so we need to collect
    // each trace (curve) number for each point by sample.
    let currentCurveNumber = traces.length;
    let curveNumbersBySample = {};
    let curveAxisBySample = {};
    let scatterDataByMetric = [];
    metrics.map((metric, metricIdx) => {
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;
      let outliers = outliersByMetric[metric];
      let valuesBySample = valuesBySampleByMetric[metric];

      let scatterData = [];
      Object.entries(valuesBySample).map(([sample, value]) => {
        if (outliers !== undefined && !outliers.includes(sample)) return; // showing only outliers

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

    let seed = 1;
    function random() {
      // Math.random does not have a seed, so we use this
      let x = Math.sin(seed++) * 10000;
      return x - Math.floor(x);
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
        let params = JSON.parse(JSON.stringify(scatterTraceParams)); // deep copy

        let color = "black"; // trace_params["marker"]["color"];
        let size = params.marker.size;
        if (highlightingEnabled) {
          color = state.highlight ?? "#cccccc";
          size = state.highlight !== null ? 10 : size;
        }

        let customData = {
          curveNumbers: curveNumbersBySample[sample],
          curveAxis: curveAxisBySample[sample],
        };

        const jitter = 0.3;
        scatters.push({
          type: "scatter",
          x: [value],
          y: [metricIdx + random() * jitter - jitter / 2], // add vertical jitter
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

  afterPlotCreated() {
    let target = this.target;
    let plot = document.getElementById(target);

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
        Plotly.Fx.hover(target, points, curveAxis);
      }
    });
  }
}
