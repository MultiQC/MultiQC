// Golden RGB contract, must stay identical to its twins in
// multiqc/plots/seqcontent.py::bin_rgb and
// templates/echarts/src/js/plots/seqcontent.js:
// R = %T, G = %A, B = %C (%G implied by the complement of the other three).
function seqContentBinRgb(bin) {
  return [Math.round((bin.t / 100) * 255), Math.round((bin.a / 100) * 255), Math.round((bin.c / 100) * 255)];
}

class SeqContentPlot extends Plot {
  constructor(dump) {
    super(dump);
    // Anchor of the auxiliary per-sample line plot built by
    // SeqContentPlot._build_drilldown_plot() (multiqc/plots/seqcontent.py), registered
    // in mqc_plots like any other plot. Absent (null) for flat reports: add_to_report()
    // skips the drilldown block entirely when self.flat.
    this.drilldownAnchor = dump["drilldown_anchor"] ?? null;
  }

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
    // Stashed for the drilldown (row index -> sample, Prev/Next): kept in lockstep with
    // filtSampleSettings, recomputed on every render so toolbox hide/rename/highlight
    // changes are reflected immediately, even while the drilldown is open.
    this._ddRows = rows;

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
    // Positions never covered by any bin for a sample stay fully transparent (alpha
    // 0), same as the old canvas simply not painting past a sample's last bin, but
    // without leaving a black block over the report background.
    let z = new Array(nSamples);
    let customdata = new Array(nSamples);
    for (let i = 0; i < nSamples; i++) {
      let sampleName = sampleSettings[i].name;
      let zRow = new Array(maxBp);
      let cdRow = new Array(maxBp);
      for (let pos = 0; pos < maxBp; pos++) {
        zRow[pos] = [0, 0, 0, 0];
        cdRow[pos] = [0, 0, 0, 0, `${pos + 1} bp`, sampleName];
      }
      for (let bin of rows[i]) {
        let [r, g, b] = seqContentBinRgb(bin);
        let label = bin.start === bin.end ? `${bin.start} bp` : `${bin.start}-${bin.end} bp`;
        for (let pos = bin.start; pos <= bin.end; pos++) {
          zRow[pos - 1] = [r, g, b, 255];
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

    this._updateYTicks(sampleSettings);

    let trace = {
      ...(dataset["trace_params"] ?? {}),
      type: "image",
      z: z,
      x0: 1,
      dx: 1,
      colormodel: "rgba",
      customdata: customdata,
      // Colored square emoji swatches map each %base line back to its heatmap
      // channel: T=red, C=blue, A=green, G=dark/gray (bin_rgb: R=%T, G=%A, B=%C).
      // Plotly hovertemplate can't render arbitrary HTML/CSS color, so emoji is the
      // portable stand-in.
      hovertemplate:
        "<b>%{customdata[5]}</b><br>%{customdata[4]}<br>" +
        "🟥 %T: %{customdata[0]}%<br>🟦 %C: %{customdata[1]}%<br>🟩 %A: %{customdata[2]}%<br>⬛ %G: %{customdata[3]}%" +
        "<extra></extra>",
    };
    return [trace];
  }

  // Recompute y tick density (array vs. auto, which sample names show) from the
  // CURRENT this.layout.height, and store it on this.layout.yaxis. Called both from
  // buildTraces() (initial render / toolbox re-render) and resize() (height-drag
  // handle), so more sample names appear as soon as the plot gets taller instead of
  // only after some other repaint recomputes it (T3 tick-repaint fix).
  _updateYTicks(sampleSettings) {
    let nSamples = sampleSettings.length;
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
    this.bindDrilldownControls();
    const gd = document.getElementById(this.anchor);
    if (gd) {
      if (!this._clickListenerBound) {
        this._clickListenerBound = true;
        gd.on("plotly_click", (data) => {
          if (!data.points || data.points.length === 0) return;
          // go.Image traces default to y0=0, dy=1: the clicked point's y is the row
          // index (0-based) into the SAME toolbox-filtered rows buildTraces() drew.
          this.openDrilldown(Math.round(data.points[0].y));
        });
      }
      if (!this._relayoutListenerBound) {
        this._relayoutListenerBound = true;
        // Box-drag zoom and double-click-to-reset both fire plotly_relayout with new
        // axis ranges. Recompute scaleratio for whatever range is on screen now, so
        // the zoomed-in view (or the reset full view) stretches to fill instead of
        // the aspect lock leaving a letterboxed strip. fixAspectRatio() is a no-op
        // once the ratio already matches (see its threshold check), so this can't
        // loop against the Plotly.relayout call it makes itself.
        gd.on("plotly_relayout", () => this.fixAspectRatio());
      }
    }
  }

  resize(newHeight, newWidth) {
    this.layout.height = newHeight;
    // Tick density is a function of the live pixel height, not computed once: redo it
    // here so dragging the plot taller shows more sample names immediately, using the
    // last render's filtered sample list (buildTraces() always populates this first).
    if (this.filtSampleSettings && this.filtSampleSettings.length > 0) {
      this._updateYTicks(this.filtSampleSettings);
    }
    super.resize(newHeight, newWidth);
    setTimeout(() => this.fixAspectRatio(), 0);
  }

  // ---------------------------------------------------------------------------------
  // Click-to-drilldown (T3.1, BUILD_PLAN.md section 1.6, option B). Shared logic: this
  // default-template class's methods are reused as-is by the echarts subclass
  // (EchartsSeqContentPlot extends window.SeqContentPlot), which only overrides
  // afterPlotCreated() to bind its own chart.on("click") instead of plotly_click.
  // ---------------------------------------------------------------------------------

  // jQuery wrapper for the standard `.mqc_hcplot_plotgroup` div that interactive_plot()
  // wraps every plot in (control panel + plot + "Created with MultiQC" footer): hiding/
  // showing this whole group is the "overview" side of the show/hide toggle.
  _overviewGroupDiv() {
    return $("#" + this.anchor).closest(".mqc_hcplot_plotgroup");
  }

  // Map a clicked/current row index (into the toolbox-filtered sample list) to that
  // sample's bins, rewrite the aux line plot's pairs and title, and render it.
  openDrilldown(rowIdx) {
    if (!this.drilldownAnchor) return; // flat report: no drilldown wrapper in the HTML
    if (!this.filtSampleSettings || rowIdx < 0 || rowIdx >= this.filtSampleSettings.length) return;
    this._ddRowIdx = rowIdx;
    this._renderDrilldown();
  }

  stepDrilldown(direction) {
    if (!this.filtSampleSettings || this.filtSampleSettings.length === 0) return;
    let n = this.filtSampleSettings.length;
    this._ddRowIdx = ((((this._ddRowIdx ?? 0) + direction) % n) + n) % n;
    this._renderDrilldown();
  }

  closeDrilldown() {
    $("#" + this.anchor + "_drilldown_wrapper").hide();
    this._overviewGroupDiv().show();
  }

  _renderDrilldown() {
    let ddPlot = mqc_plots[this.drilldownAnchor];
    if (!ddPlot) return;
    let sampleSetting = this.filtSampleSettings[this._ddRowIdx];
    let bins = [...this._ddRows[this._ddRowIdx]].sort((a, b) => a.start - b.start);

    // Series order fixed at creation time (multiqc/plots/seqcontent.py::_DRILLDOWN_SERIES):
    // % T, % C, % A, % G.
    const baseKeys = ["t", "c", "a", "g"];
    let dataset = ddPlot.datasets[0];
    dataset.lines.forEach((line, i) => {
      line.pairs = bins.map((b) => [b.start, b[baseKeys[i]]]);
    });

    let title = sampleSetting.name;
    // The default (Plotly) template's renderPlot() does
    // `updateObject(plot.layout, dataset.layout, false)` on every render, which
    // unconditionally overwrites plot.layout.title with dataset.layout.title; both have
    // to be updated or the dataset-level (stale) title wins.
    ddPlot.layout.title = ddPlot.layout.title || {};
    ddPlot.layout.title.text = title;
    dataset.layout.title = dataset.layout.title || {};
    dataset.layout.title.text = title;
    // ECharts skeleton (present when config.plotting_engine == "echarts"): each
    // dataset's option is a self-contained clone, no plot-level merge to fight.
    if (ddPlot.echarts) {
      let ddLayout = ddPlot.echarts.datasets[0].layout;
      ddLayout.title = ddLayout.title || {};
      ddLayout.title.text = title;
    }

    this._overviewGroupDiv().hide();
    $("#" + this.anchor + "_drilldown_wrapper").show();
    window.renderPlot(this.drilldownAnchor);
  }

  bindDrilldownControls() {
    if (this._drilldownBound || !this.drilldownAnchor) return;
    this._drilldownBound = true;
    $("#" + this.anchor + "_drilldown_back").on("click", () => this.closeDrilldown());
    $("#" + this.anchor + "_drilldown_prev").on("click", () => this.stepDrilldown(-1));
    $("#" + this.anchor + "_drilldown_next").on("click", () => this.stepDrilldown(1));
  }
}

// Make SeqContentPlot globally available
window.SeqContentPlot = SeqContentPlot;
