class ViolinPlot extends Plot {
  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.active_dataset_idx]["all_samples"].length;
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let dataset = this.datasets[this.active_dataset_idx];
    let valuesByMetric = dataset["values_by_metric"];
    let samplesByMetric = dataset["samples_by_metric"];
    if (Object.keys(valuesByMetric).length === 0) return [];
    let headersByMetric = dataset["headers_by_metric"];
    let outlierIndicesByMetric = dataset["outlier_indices_by_metric"];
    let layout = this.layout;
    const trace_params = this.trace_params;

    let allSamples = this.datasets[this.active_dataset_idx]["all_samples"];
    let sampleSettings = applyToolboxSettings(allSamples);
    if (sampleSettings == null) return; // All series are hidden, do not render the graph.

    let filteredValuesByMetric = {};
    let filteredSamplesByMetric = {};
    Object.entries(valuesByMetric).map(([metric, values]) => {
      filteredValuesByMetric[metric] = [];
      filteredSamplesByMetric[metric] = [];
      let samples = samplesByMetric[metric];
      for (let i = 0; i < samples.length; i++) {
        let sample = samples[i];
        if (sampleSettings[samples.indexOf(sample)].hidden) continue;
        filteredValuesByMetric[metric].push(values[i]);
        filteredSamplesByMetric[metric].push(sample);
      }
    });
    valuesByMetric = filteredValuesByMetric;
    samplesByMetric = filteredSamplesByMetric;

    let traces = [];
    Object.keys(valuesByMetric).map((metric, metricIdx) => {
      let params = JSON.parse(JSON.stringify(trace_params)); // deep copy
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;

      let values = valuesByMetric[metric];
      let samples = samplesByMetric[metric];
      let header = headersByMetric[metric];

      traces.push({
        type: "violin",
        x: values,
        name: metricIdx, // headerByMetric[metric].title + "  ",
        text: samples, // sample names
        xaxis: "x" + axisKey,
        yaxis: "y" + axisKey,
        ...params,
      });

      let title = header.title + "  ";
      if (header["namespace"] !== null) title = header["namespace"] + "  <br>" + title;

      // Set layouts for each violin individually
      layout["xaxis" + axisKey] = Object.assign(JSON.parse(JSON.stringify(layout.xaxis)), header["xaxis"]);
      layout["yaxis" + axisKey] = Object.assign(JSON.parse(JSON.stringify(layout.yaxis)), {
        tickmode: "array",
        tickvals: [metricIdx],
        ticktext: [title],
      });

      if (header["color"] !== null) {
        let color = "rgba(" + header["color"] + ",1)";
        layout["yaxis" + axisKey]["tickfont"] = { color: color };
      }
    });

    // We want to select points on all violins belonging to this specific sample.
    // Points are rendered each as a separate trace, so we need to collect
    // each trace (curve) number for each point by sample.
    // TODO: a potentially more efficient alternative, try to add all points
    //  for a sample as a single trace?
    let currentCurveNumber = Object.keys(valuesByMetric).length;
    let curveNumbersBySample = {};
    let curveAxisBySample = {};
    let scatterDataByMetric = [];
    Object.keys(valuesByMetric).map((metric, metricIdx) => {
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;
      let values = valuesByMetric[metric];
      let samples = samplesByMetric[metric];
      let outlierIndices = outlierIndicesByMetric[metric];

      let scatterValues = [];
      for (let sampleIdx = 0; sampleIdx < samples.length; sampleIdx++) {
        if (outlierIndices !== undefined && !outlierIndices.includes(sampleIdx)) continue; // skip non-outliers

        let sample = samples[sampleIdx];
        let value = values[sampleIdx];
        scatterValues.push([sample, value]);

        if (!curveNumbersBySample[sample]) {
          curveNumbersBySample[sample] = [];
          curveAxisBySample[sample] = [];
        }
        curveNumbersBySample[sample].push(currentCurveNumber++);
        curveAxisBySample[sample].push("x" + axisKey + "y" + axisKey);
      }
      scatterDataByMetric.push(scatterValues);
    });

    let highlightingEnabled = sampleSettings.filter((s) => s.highlight).length > 0;

    let scatters = [];
    // Add scatter plots on top of violins to show individual points
    // Plotly supports showing points automatically with `points="all"`,
    // however, it's problematic to give each sample individual color,
    // and set up hover events to highlight all points for a sample. So as
    // a workaround, we add a separate scatter plot. One thing to solve later
    // would be to be able to only show outliers, not all points, as the violin
    // plot can do that automatically, and we lose this functionality here.
    scatterDataByMetric.map((sampleValues, metricIdx) => {
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;

      sampleValues.map(([sample, value]) => {
        let state = sampleSettings[allSamples.indexOf(sample)];
        let params = JSON.parse(JSON.stringify(trace_params)); // deep copy

        let color = "black"; // trace_params["marker"]["color"];
        if (highlightingEnabled) color = state.highlight ?? "#cccccc";

        let customData = {
          curveNumbers: curveNumbersBySample[sample],
          curveAxis: curveAxisBySample[sample],
          defaultMarker: {
            color: color,
            size: state.highlight !== null ? 10 : 4,
            line: {},
          },
        };

        const jitter = 0.5;
        scatters.push({
          type: "scatter",
          mode: "markers",
          x: [value],
          y: [metricIdx + Math.random() * jitter - jitter / 2], // add vertical jitter
          text: [state.name ?? sample],
          xaxis: "x" + axisKey,
          yaxis: "y" + axisKey,
          showlegend: false,
          hovertemplate: params["hovertemplate"],
          hoverlabel: { bgcolor: "white" },
          customdata: customData,
          marker: JSON.parse(JSON.stringify(customData.defaultMarker)),
          ...params,
        });
      });
    });
    traces = traces.concat(scatters);
    return traces;
  }

  // getMarkerStyle(highlighted) {
  //   return {
  //     color: color,
  //     size: settings.highlight !== null ? 10 : 4,
  //     line: { width: 0 },
  //   };
  // }

  afterPlotCreated() {
    let target = this.target;
    let trace_params = this.trace_params; // deep copy
    let plot = document.getElementById(target);

    plot
      .on("plotly_hover", function (eventdata) {
        console.log("hover");
        if (!eventdata.points) return;
        let point = eventdata.points[0];
        // if (point.data.type === "violin") {
        //   console.log("hover violin", point);
        //
        //   let curveNumbers = point.data.customdata["curveNumbers"][point.pointNumber];
        //   let curveAxes = point.data.customdata["curveAxis"][point.pointNumber];
        //   let points = curveNumbers.map((curveNum) => {
        //     return { curveNumber: curveNum, pointNumber: 0 };
        //   });
        //   Plotly.Fx.hover(target, points, curveAxes);
        // }

        if (point.data.type === "scatter") {
          let update = {
            "marker.size": 10,
            "marker.color": "rgb(98,192,248)",
            "marker.line.color": "black",
            "marker.line.width": 2,
          };
          Plotly.restyle(target, update, point.data.customdata["curveNumbers"]);

          // let curveNumbers = point.data.customdata["curveNumbers"];
          // let curveAxis = point.data.customdata["curveAxis"];
          // let points = curveNumbers.map((curveNum) => {
          //   return { curveNumber: curveNum, pointNumber: 0 };
          // });
          // Plotly.Fx.hover(target, points, curveAxis);
        }
      })
      .on("plotly_unhover", function (eventdata) {
        console.log("unhover");
        if (!eventdata.points) return;
        let point = eventdata.points[0];
        if (point.data.type === "scatter") {
          let curveNumbers = point.data.customdata["curveNumbers"];
          let marker = JSON.parse(JSON.stringify(point.data.customdata["defaultMarker"]));
          let update = { marker: marker };
          Plotly.restyle(target, update, curveNumbers);
        }
      });
  }
}
