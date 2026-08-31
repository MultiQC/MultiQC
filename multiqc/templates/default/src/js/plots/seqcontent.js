// Golden RGB contract, must stay identical to its twins in
// multiqc/plots/seqcontent.py::bin_rgb and
// templates/echarts/src/js/plots/seqcontent.js:
// R = %T, G = %A, B = %C (%G implied by the complement of the other three).
function seqContentBinRgb(bin) {
  return [Math.round((bin.t / 100) * 255), Math.round((bin.a / 100) * 255), Math.round((bin.c / 100) * 255)];
}

class SeqContentPlot extends Plot {
  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    return this.datasets[this.activeDatasetIdx].samples.length; // no samples
  }

  prepData(dataset) {
    // Prepare data to either build Plotly traces or export as a file
    dataset = dataset ?? this.datasets[this.activeDatasetIdx];

    let sampleSettings = applyToolboxSettings(dataset.samples);
    let rows = dataset.rows.filter((row, idx) => !sampleSettings[idx].hidden);
    this.filtSampleSettings = sampleSettings.filter((s) => !s.hidden);

    return [this.filtSampleSettings, rows];
  }

  plotAiHeader() {
    let result = super.plotAiHeader();
    if (this.pconfig.xlab) result += `X axis: ${this.pconfig.xlab}\n`;
    return result;
  }

  formatDatasetForAiPrompt(dataset) {
    let [sampleSettings, rows] = this.prepData(dataset);

    if (sampleSettings.length === 0) {
      return "All samples are hidden by user, so no data to analyse. Please inform user to use the toolbox to unhide samples.\n";
    }

    let prompt = "";
    sampleSettings.forEach((s, i) => {
      let name = s.pseudonym ?? s.name;
      let positions = rows[i].map((b) => `${b.label}: T=${b.t}/C=${b.c}/A=${b.a}/G=${b.g}`).join(", ");
      prompt += `${name}: ${positions}\n`;
    });
    return prompt;
  }

  buildTraces() {
    let [sampleSettings, rows] = this.prepData();
    if (sampleSettings.length === 0) return [];

    let dataset = this.datasets[this.activeDatasetIdx];
    let maxBp = Math.max(dataset.max_bp, 1);
    let nSamples = rows.length;

    // Expand the (possibly non-uniform) bins into a uniform per-sample pixel row.
    // Positions never covered by any bin for a sample stay black, same as the old
    // canvas simply not painting past a sample's last bin.
    let z = new Array(nSamples);
    let customdata = new Array(nSamples);
    for (let i = 0; i < nSamples; i++) {
      let sampleName = sampleSettings[i].name;
      let zRow = new Array(maxBp);
      let cdRow = new Array(maxBp);
      for (let pos = 0; pos < maxBp; pos++) {
        zRow[pos] = [0, 0, 0];
        cdRow[pos] = [0, 0, 0, 0, `${pos + 1} bp`, sampleName];
      }
      for (let bin of rows[i]) {
        let [r, g, b] = seqContentBinRgb(bin);
        let label = bin.start === bin.end ? `${bin.start} bp` : `${bin.start}-${bin.end} bp`;
        for (let pos = bin.start; pos <= bin.end; pos++) {
          zRow[pos - 1] = [r, g, b];
          cdRow[pos - 1] = [bin.t, bin.c, bin.a, bin.g, label, sampleName];
        }
      }
      z[i] = zRow;
      customdata[i] = cdRow;
    }

    // y tick labels are the sample names (heatmap.js pattern); an image trace's y
    // coordinate is the row index (0..n-1), not the sample name itself, so tick
    // positions and tick labels have to be tracked separately.
    this.layout.yaxis.tickmode = "array";
    this.layout.yaxis.tickvals = [...Array(nSamples).keys()];
    this.layout.yaxis.ticktext = sampleSettings.map((s) => s.name);

    const maxYTicks = (this.layout.height - 200) / 12;
    this.recalculateTicks(sampleSettings, this.layout.yaxis, maxYTicks);
    if (this.layout.yaxis.tickvals) {
      // recalculateTicks positions ticks using the sample name itself (correct for
      // heatmap's categorical axis, where the category value IS the sample name).
      // Our axis is numeric row indices, so remap the generated tickvals from
      // names back to their row index, keeping the colored ticktext it built.
      this.layout.yaxis.tickvals = this.layout.yaxis.tickvals.map((name) =>
        sampleSettings.findIndex((s) => s.name === name),
      );
    } else {
      // No highlighted samples: recalculateTicks reverts to Plotly's "auto" ticks,
      // which would show plain numeric row indices on this axis instead of sample
      // names. Categorical axes (heatmap.js) get free tick-skipping from Plotly
      // when auto; our axis is numeric array-mode, so thin the ticks ourselves to
      // roughly maxYTicks, evenly spaced, to avoid every sample label overlapping.
      let step = Math.max(1, Math.ceil(nSamples / Math.max(1, maxYTicks)));
      let idxs = [];
      for (let i = 0; i < nSamples; i += step) idxs.push(i);
      if (idxs[idxs.length - 1] !== nSamples - 1) idxs.push(nSamples - 1);

      this.layout.yaxis.tickmode = "array";
      this.layout.yaxis.tickvals = idxs;
      this.layout.yaxis.ticktext = idxs.map((i) => sampleSettings[i].name);
    }

    let trace = {
      ...(dataset["trace_params"] ?? {}),
      type: "image",
      z: z,
      x0: 1,
      dx: 1,
      colormodel: "rgb",
      customdata: customdata,
      hovertemplate:
        "<b>%{customdata[5]}</b><br>%{customdata[4]}<br>" +
        "%T: %{customdata[0]}%<br>%C: %{customdata[1]}%<br>%A: %{customdata[2]}%<br>%G: %{customdata[3]}%" +
        "<extra></extra>",
    };
    return [trace];
  }

  exportData(format) {
    let [sampleSettings, rows] = this.prepData();
    let sep = format === "tsv" ? "\t" : ",";

    let csv = ["Sample", "Position", "A", "C", "G", "T"].join(sep) + "\n";
    sampleSettings.forEach((s, i) => {
      rows[i].forEach((bin) => {
        csv += [s.name, bin.label, bin.a, bin.c, bin.g, bin.t].join(sep) + "\n";
      });
    });
    return csv;
  }

  // Image traces lock the y axis to an equal-aspect "1 data unit == 1 screen px in
  // both axes" via yaxis.scaleanchor, and that coupling cannot be cleared (verified
  // against plotly.js 3.1.2: setting scaleanchor to null/false in the layout, at
  // creation or via a later relayout, is silently overwritten back to "x"). The
  // supported knob to distort that forced aspect is yaxis.scaleratio, which we set
  // to exactly the ratio that makes both axis domains fill the whole plot area, so
  // the image stretches to fill like the old canvas instead of leaving a
  // letterboxed strip.
  fixAspectRatio() {
    const gd = document.getElementById(this.anchor);
    if (!gd || !gd._fullLayout) return;
    const fl = gd._fullLayout;
    if (!fl.xaxis || !fl.yaxis || !fl._size || !fl._size.w || !fl._size.h) return;

    const xRange = Math.abs(fl.xaxis.range[1] - fl.xaxis.range[0]);
    const yRange = Math.abs(fl.yaxis.range[1] - fl.yaxis.range[0]);
    if (!xRange || !yRange) return;

    const scaleratio = (fl._size.h / fl._size.w) * (xRange / yRange);
    if (Math.abs((fl.yaxis.scaleratio ?? 1) - scaleratio) < 0.01) return;

    this.layout.yaxis.scaleratio = scaleratio;
    Plotly.relayout(this.anchor, { "yaxis.scaleratio": scaleratio });
  }

  afterPlotCreated() {
    this.fixAspectRatio();
    if (!this._resizeListenerBound) {
      this._resizeListenerBound = true;
      window.addEventListener("resize", () => {
        clearTimeout(this._resizeTimer);
        this._resizeTimer = setTimeout(() => this.fixAspectRatio(), 100);
      });
    }
  }

  resize(newHeight, newWidth) {
    super.resize(newHeight, newWidth);
    setTimeout(() => this.fixAspectRatio(), 0);
  }
}

// Make SeqContentPlot globally available
window.SeqContentPlot = SeqContentPlot;
