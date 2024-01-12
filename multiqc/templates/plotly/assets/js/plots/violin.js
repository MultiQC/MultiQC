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

    let violins = [];
    Object.keys(valuesByMetric).map((metric, metricIdx) => {
      let params = JSON.parse(JSON.stringify(trace_params)); // deep copy
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;

      let values = valuesByMetric[metric];
      let samples = samplesByMetric[metric];
      let header = headersByMetric[metric];

      if (header["color"] !== null) params["fillcolor"] = "rgba(" + header["color"] + ",0.5)";

      violins.push({
        type: "violin",
        x: values,
        name: metricIdx, // headerByMetric[metric].title + "  ",
        text: samples, // sample names
        xaxis: "x" + axisKey,
        yaxis: "y" + axisKey,
        ...params,
      });

      // Set layouts for each violin individually
      layout["xaxis" + axisKey] = Object.assign(JSON.parse(JSON.stringify(layout.xaxis)), header["xaxis"]);
      layout["yaxis" + axisKey] = Object.assign(JSON.parse(JSON.stringify(layout.yaxis)), {
        tickmode: "array",
        tickvals: [metricIdx],
        ticktext: [header.title + "  "],
      });
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

    let highlighting = sampleSettings.filter((s) => s.highlight).length > 0;

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
        let settings = sampleSettings[allSamples.indexOf(sample)];
        let params = JSON.parse(JSON.stringify(trace_params)); // deep copy

        let color = trace_params["marker"]["color"];
        if (highlighting) color = settings.highlight ?? "#cccccc";

        const jitter = 0.5;
        scatters.push({
          type: "scatter",
          mode: "markers",
          x: [value],
          y: [metricIdx + Math.random() * jitter - jitter / 2], // add vertical jitter
          text: [settings.name ?? sample],
          xaxis: "x" + axisKey,
          yaxis: "y" + axisKey,
          highlighted: settings.highlight !== null,
          marker: { color: color, size: settings.highlight !== null ? 10 : 6 },
          showlegend: false,
          hovertemplate: params["hovertemplate"],
          customdata: { curveNumbers: curveNumbersBySample[sample], curveAxis: curveAxisBySample[sample] },
          ...params,
        });
      });
    });
    return violins.concat(scatters);
  }

  afterPlotCreated() {
    let target = this.target;
    let trace_params = this.trace_params; // deep copy
    let plot = document.getElementById(target);

    plot.on("plotly_hover", function (eventdata) {
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
        console.log("hover scatter", point);

        // let update = {
        //   "marker.size": 10,
        //   "marker.color": "rgb(55,126,184)",
        //   "marker.line.color": "black",
        //   "marker.line.width": 1,
        // };

        let curveNumbers = point.data.customdata["curveNumbers"];
        let curveAxis = point.data.customdata["curveAxis"];
        let points = curveNumbers.map((curveNum) => {
          return { curveNumber: curveNum, pointNumber: 0 };
        });
        Plotly.Fx.hover(target, points, curveAxis);
      }
    });
    // .on("plotly_unhover", function (eventdata) {
    //   if (!eventdata.points) return;
    //   let point = eventdata.points[0];
    //   if (point.data.type === "scatter") {
    //     let marker = JSON.parse(JSON.stringify(trace_params["marker"]));
    //     let curveNumbers = point.data.customdata["curveNumbers"];
    //     console.log("unhover", point);
    //     // TODO: only change the marker size back if the point is not highlighted
    //     let update = {
    //       marker: marker,
    //     };
    //     Plotly.restyle(target, update, curveNumbers);
    //   }
    // });
  }
}
