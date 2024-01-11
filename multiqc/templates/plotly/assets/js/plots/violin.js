class ViolinPlot extends Plot {
  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    if (this.datasets[this.active_dataset_idx]["samples"].length === 0) return 0; // no samples in a dataset
    return this.datasets[this.active_dataset_idx]["header_by_metric"].length; // no metrics in a dataset
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let dataByMetric = this.datasets[this.active_dataset_idx]["data_by_metric"];
    if (dataByMetric.length === 0) return [];
    let headerByMetric = this.datasets[this.active_dataset_idx]["header_by_metric"];
    let layout = this.layout;

    let samples = this.datasets[this.active_dataset_idx]["samples"];
    let sampleSettings = applyToolboxSettings(samples);
    if (sampleSettings == null) return; // All series are hidden, do not render the graph.

    let filteredDataByMetric = [];
    Object.entries(dataByMetric).map(([metric, valueBySample]) => {
      filteredDataByMetric[metric] = {};
      Object.entries(valueBySample).map(([sample, value]) => {
        if (!sampleSettings[samples.indexOf(sample)].hidden) {
          filteredDataByMetric[metric][sample] = value;
        }
      });
    });
    dataByMetric = filteredDataByMetric;

    let violins = [];
    Object.keys(dataByMetric).map((metric, metricIdx) => {
      let params = JSON.parse(JSON.stringify(this.trace_params)); // deep copy
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;

      // Set layouts for each violin individually
      layout["xaxis" + axisKey] = Object.assign(
        JSON.parse(JSON.stringify(layout.xaxis)),
        headerByMetric[metric]["xaxis"],
      );
      layout["yaxis" + axisKey] = JSON.parse(JSON.stringify(layout.yaxis));

      let valueBySample = dataByMetric[metric];
      violins.push({
        type: "violin",
        x: Object.values(valueBySample),
        name: headerByMetric[metric].title + "  ",
        text: Object.keys(valueBySample), // sample names
        xaxis: "x" + axisKey,
        yaxis: "y" + axisKey,
        ...params,
      });
    });

    // We want to select points on all violins belonging to this specific sample.
    // Points are rendered each as a separate trace, so we need to collect
    // each trace (curve) number for each point by sample.
    // TODO: a potentially more efficient alternative, try to add all points
    //  for a sample as a single trace?
    let currentCurveNumber = violins.length;
    let curveNumbersBySample = {};
    let curveAxisBySample = {};
    Object.keys(dataByMetric).map((metric, metricIdx) => {
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;
      let valueBySample = dataByMetric[metric];
      Object.keys(valueBySample).map((sample) => {
        if (!curveNumbersBySample[sample]) {
          curveNumbersBySample[sample] = [];
          curveAxisBySample[sample] = [];
        }
        curveNumbersBySample[sample].push(currentCurveNumber++);
        curveAxisBySample[sample].push("x" + axisKey + "y" + axisKey);
      });
    });

    // Add scatter plots on top of violins to show individual points
    // Plotly supports showing points automatically with `points="all"`,
    // however, it's problematic to give each sample individual color,
    // and set up hover events to highlight all points for a sample. So as
    // a workaround, we add a separate scatter plot. One thing to solve later
    // would be to be able to only show outliers, not all points, as the violin
    // plot can do that automatically, and we lose this functionality here.
    let points = [];

    let highlighting = sampleSettings.filter((s) => s.highlight).length > 0;

    Object.keys(dataByMetric).map((metric, idx) => {
      let params = JSON.parse(JSON.stringify(this.trace_params)); // deep copy
      let axisKey = idx === 0 ? "" : idx + 1;

      Object.entries(dataByMetric[metric]).map(([sample, value]) => {
        let sampleData = sampleSettings[samples.indexOf(sample)];

        let color = "rgb(55,126,184)";
        if (highlighting) color = sampleData.highlight ?? "#cccccc";
        if (sampleData.highlight) {
          console.log("highlighted");
        }

        let sampleTrace = {
          type: "scatter",
          mode: "markers",
          x: [value],
          y: [headerByMetric[metric].title + "  "],
          text: [sampleData.name ?? sample],
          xaxis: "x" + axisKey,
          yaxis: "y" + axisKey,
          highlighted: sampleData.highlight !== null,
          marker: { color: color, size: sampleData.highlight !== null ? 10 : 6 },
          showlegend: false,
          hovertemplate: params["hovertemplate"],
          customdata: { curveNumbers: curveNumbersBySample[sample], curveAxis: curveAxisBySample[sample] },
          ...params,
        };
        points.push(sampleTrace);
      });
    });

    // Reorder points so highlighted points are on top
    // let highlighted = points.filter((trace) => trace.highlighted);
    // let nonHighlighted = points.filter((trace) => !trace.highlighted);
    // points = nonHighlighted.concat(highlighted);

    return violins.concat(points);
  }

  afterPlotCreated(target) {
    let plot = document.getElementById(target);

    plot
      .on("plotly_hover", function (eventdata) {
        if (!eventdata.points) return;
        let point = eventdata.points[0];
        if (point.data.type === "scatter") {
          let curveNumbers = point.data.customdata["curveNumbers"];
          // let curveAxis = point.data.customdata["curveAxis"];
          // let points = curveNumbers.map((curveNum) => {
          //   return { curveNumber: curveNum, pointNumber: 0 };
          // });
          // Plotly.Fx.hover(target, points, curveAxis);
          let update = {
            "marker.size": 10,
            "marker.line.color": "black",
            "marker.line.width": 1,
          };
          Plotly.restyle(target, update, curveNumbers);
        }
      })
      .on("plotly_unhover", function (eventdata) {
        // TODO: only change the marker size back if the point is not highlighted
        let color = "rgb(55,126,184)";
        Plotly.restyle(target, {
          "marker.size": 6,
          "marker.line.color": color,
          "marker.line.width": 0,
        });
      });
  }
}
