class ViolinPlot extends Plot {
  constructor(dump) {
    super(dump);
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.active_dataset_idx]["data_by_metric"].length; // no samples in a dataset
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let dataByMetric = this.datasets[this.active_dataset_idx]["data_by_metric"];
    if (dataByMetric.length === 0) return [];
    let headerByMetric = this.datasets[this.active_dataset_idx]["header_by_metric"];
    let sampleColors = this.datasets[this.active_dataset_idx]["sample_colors"];
    let layout = this.layout;

    let violins = Object.keys(headerByMetric).map((metric, idx) => {
      let params = JSON.parse(JSON.stringify(this.trace_params)); // deep copy
      let axisKey = idx === 0 ? "" : idx + 1;

      // set layouts for each violin
      layout["xaxis" + axisKey] = Object.assign(
        JSON.parse(JSON.stringify(layout.xaxis)),
        headerByMetric[metric]["xaxis"],
      );
      layout["yaxis" + axisKey] = JSON.parse(JSON.stringify(layout.yaxis));

      return {
        type: "violin",
        x: Object.values(dataByMetric[metric]),
        name: headerByMetric[metric].title + "  ",
        text: Object.keys(dataByMetric[metric]), // sample names
        xaxis: "x" + axisKey,
        yaxis: "y" + axisKey,
        ...params,
      };
    });

    // Add scatter plots on top of violins to show individual points
    // Plotly supports showing points automatically with `points="all"`,
    // however, it's problematic to give each sample individual color,
    // and set up hover events to highlight all points for a sample. So as
    // a workaround, we add a separate scatter plot. One thing to solve later
    // would be to be able to only show outliers, not all points, as the violin
    // plot can do that automatically, and we lose this functionality here.
    // ---
    // We want to select points on all violins belonging to this specific sample.
    // Points are rendered each as a separate trace, so we need to collect
    // each trace (curve) number for each point by sample.
    // TODO: a potentially more efficient alternative, try to add all points
    //  for a sample as a single trace?
    let currentCurveNumber = violins.length;
    let curveNumbersBySample = {};
    let curveAxisBySample = {};
    Object.keys(headerByMetric).map((metric, metricIdx) => {
      let axisKey = metricIdx === 0 ? "" : metricIdx + 1;
      Object.entries(dataByMetric[metric]).map(([sample, value]) => {
        if (!curveNumbersBySample[sample]) {
          curveNumbersBySample[sample] = [];
          curveAxisBySample[sample] = [];
        }
        curveNumbersBySample[sample].push(currentCurveNumber++);
        curveAxisBySample[sample].push("x" + axisKey + "y" + axisKey);
      });
    });

    let points = [];
    Object.keys(headerByMetric).map((metric, idx) => {
      let params = JSON.parse(JSON.stringify(this.trace_params)); // deep copy
      let axisKey = idx === 0 ? "" : idx + 1;

      Object.entries(dataByMetric[metric]).map(([sample, value]) => {
        let sampleColor = sampleColors[sample];
        let sampleTrace = {
          type: "scatter",
          mode: "markers",
          x: [value],
          y: [headerByMetric[metric].title + "  "],
          text: [sample],
          xaxis: "x" + axisKey,
          yaxis: "y" + axisKey,
          marker: {
            color: sampleColor,
          },
          showlegend: false,
          customdata: { curveNumbers: curveNumbersBySample[sample], curveAxis: curveAxisBySample[sample] },
          ...params,
        };
        points.push(sampleTrace);
      });
    });
    return violins.concat(points);
  }

  afterPlotCreated(target) {
    let plot = document.getElementById(target);
    plot.on("plotly_hover", function (eventdata) {
      let point = eventdata.points[0];
      if (point.data.type === "scatter") {
        let curveNumbers = point.data.customdata["curveNumbers"];
        let curveAxis = point.data.customdata["curveAxis"];
        console.log("plotly_click", point, curveNumbers, curveAxis);
        console.log("plotly_click: triggering hover of points", curveNumbers, curveAxis);
        let points = curveNumbers.map((curveNum) => {
          return { curveNumber: curveNum, pointNumber: 0 };
        });
        Plotly.Fx.hover(target, points, curveAxis);
      }
      console.log("plotly_click: finished");
    });
  }
}
