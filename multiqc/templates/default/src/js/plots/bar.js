class BarPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.filteredSettings = [];
  }

  activeDatasetSize() {
    if (this.datasets.length === 0) return 0; // no datasets
    let cats = this.datasets[this.activeDatasetIdx]["cats"];
    if (cats.length === 0) return 0; // no categories
    return cats[0].data.length; // no data for a category
  }

  prepData(dataset) {
    dataset = dataset ?? this.datasets[this.activeDatasetIdx];
    let cats = dataset["cats"];
    let samples = dataset["samples"];
    let groupLabels = dataset["group_labels"] || null;
    let offsetGroups = dataset["offset_groups"] || null;

    let samplesSettings = applyToolboxSettings(samples);

    // Rename and filter samples:
    this.filteredSettings = samplesSettings.filter((s) => !s.hidden);

    // Filter group labels to match visible samples
    let filteredGroupLabels = null;
    if (groupLabels) {
      filteredGroupLabels = groupLabels.filter((_, si) => !samplesSettings[si].hidden);
    }
    this.filteredGroupLabels = filteredGroupLabels;
    this.offsetGroups = offsetGroups;

    cats = cats.map((cat) => {
      let data = this.pActive ? cat["data_pct"] : cat.data;
      return {
        data: data.filter((_, si) => !samplesSettings[si].hidden),
        color: cat.color, // formatted as "r,g,b", to be wrapped with "rgb()" or "rgba()"
        name: cat.name,
      };
    });

    return [cats];
  }

  plotAiHeader() {
    let result = super.plotAiHeader();
    if (this.pconfig.ylab) result += `Values: ${this.pconfig.ylab}\n`;
    return result;
  }

  formatDatasetForAiPrompt(dataset) {
    let prompt = "";

    let cats = dataset.cats;
    let samples = dataset.samples;
    let samplesSettings = applyToolboxSettings(samples);

    // Check if all samples are hidden
    if (samplesSettings.every((s) => s.hidden)) {
      prompt +=
        "All samples are hidden by user, so no data to analyse. Please inform user to use the toolbox to unhide samples.\n";
      return prompt;
    }

    prompt += "|Sample|" + cats.map((cat) => cat.name).join("|") + "|\n";
    prompt += "|---|" + cats.map(() => "---").join("|") + "|\n";

    let suffix = "";
    if (this.pActive) {
      suffix += "%";
      if (this.layout.xaxis.ticksuffix && this.layout.xaxis.ticksuffix !== "%") {
        suffix += " " + this.layout.xaxis.ticksuffix;
      }
    } else if (this.layout.xaxis.ticksuffix) {
      suffix += " " + this.layout.xaxis.ticksuffix;
    }

    // Create data rows

    samplesSettings.forEach((sample, idx) => {
      if (sample.hidden) return;
      prompt +=
        `|${sample.pseudonym ?? sample.name}|` +
        cats
          .map((cat) => {
            let val = this.pActive ? cat.data_pct[idx] : cat.data[idx];
            val = !Number.isFinite(val) ? "" : Number.isInteger(val) ? val : parseFloat(val.toFixed(2));
            if (val !== "" && suffix) val += suffix;
            return val;
          })
          .join("|") +
        "|\n";
    });

    return prompt;
  }

  resize(newHeight) {
    this.layout.height = newHeight;

    // Only recalculate ticks for non-multicategory axes
    let useMulticategory = this.filteredGroupLabels && this.filteredGroupLabels.length > 0;
    if (!useMulticategory) {
      const maxTicks = (this.layout.height - 140) / 12;
      this.recalculateTicks(this.filteredSettings, this.layout.yaxis, maxTicks);
    }

    super.resize(newHeight);
  }

  // Constructs and returns traces for the Plotly plot
  buildTraces() {
    let [cats] = this.prepData();
    if (cats.length === 0 || this.filteredSettings.length === 0) return [];

    // Check if we have group labels for multicategory axis
    let useMulticategory = this.filteredGroupLabels && this.filteredGroupLabels.length > 0;

    // Only recalculate ticks for non-multicategory axes
    // (multicategory axes have their own tick handling that conflicts with tickvals/ticktext)
    if (!useMulticategory) {
      const maxTicks = (this.layout.height - 140) / 12;
      this.recalculateTicks(this.filteredSettings, this.layout.yaxis, maxTicks);
    }

    let highlighted = this.filteredSettings.filter((s) => s.highlight);
    let firstHighlightedSample = this.firstHighlightedSample(this.filteredSettings);
    let traceParams = this.datasets[this.activeDatasetIdx]["trace_params"];

    // When using sample groups with offsetgroup, create traces per (category, sample) combination
    // This creates natural gaps between groups using Plotly's offsetgroup feature
    // y-axis shows group labels (e.g., read lengths), offsetgroup creates sample sub-groups
    if (useMulticategory) {
      // Single pass: build all data structures we need
      let uniqueSamples = [];
      let seenSamples = new Set();
      let sampleGroupData = {}; // sampleName -> groupLabel -> [{name, value}, ...]
      let sampleHighlightStatus = {}; // sampleName -> boolean (is any entry highlighted?)
      let sampleGroupEntries = {}; // sampleName -> [{groupLabel, dataIdx}, ...]

      this.filteredSettings.forEach((sample, idx) => {
        let displayName = sample.name;
        let groupLabel = this.filteredGroupLabels[idx];

        // Track unique samples
        if (!seenSamples.has(displayName)) {
          seenSamples.add(displayName);
          uniqueSamples.push(displayName);
          sampleGroupData[displayName] = {};
          sampleHighlightStatus[displayName] = false;
          sampleGroupEntries[displayName] = [];
        }

        // Track highlight status
        if (sample.highlight) {
          sampleHighlightStatus[displayName] = true;
        }

        // Store group entry for this sample
        sampleGroupEntries[displayName].push({ groupLabel, dataIdx: idx });

        // Build category data for this sample+group combination
        if (!sampleGroupData[displayName][groupLabel]) {
          sampleGroupData[displayName][groupLabel] = cats.map((cat) => ({
            name: cat.name,
            value: cat.data[idx],
          }));
        }
      });

      // Get hover format from dataset layout, fallback to .2f
      let hoverFormat = this.layout.xaxis?.hoverformat || ".2f";

      // Build hovertemplate once per category (reusable across samples)
      let hoverTemplates = {};
      cats.forEach((cat) => {
        let hoverLines = ["<b>%{customdata.sampleName}</b>"];
        cats.forEach((c, ci) => {
          let prefix = c.name === cat.name ? "► <b>" : "   ";
          let suffix = c.name === cat.name ? "</b>" : "";
          hoverLines.push(prefix + c.name + ": %{customdata.catData[" + ci + "].value:" + hoverFormat + "}" + suffix);
        });
        hoverTemplates[cat.name] = hoverLines.join("<br>") + "<extra></extra>";
      });

      // Create traces
      let traces = [];
      uniqueSamples.forEach((sampleName, sampleIdx) => {
        // Determine alpha once per sample
        let sampleHighlighted = sampleHighlightStatus[sampleName];
        let alpha = highlighted.length > 0 && !sampleHighlighted ? 0.1 : 1;

        cats.forEach((cat) => {
          // Collect data for this sample across all its groups
          let groupLabels = [];
          let groupData = [];
          let customData = [];

          sampleGroupEntries[sampleName].forEach(({ groupLabel, dataIdx }) => {
            groupLabels.push(groupLabel);
            groupData.push(cat.data[dataIdx]);
            customData.push({
              sampleName: sampleName,
              catData: sampleGroupData[sampleName][groupLabel],
            });
          });

          if (groupLabels.length > 0) {
            traces.push({
              type: "bar",
              x: groupData,
              y: groupLabels,
              customdata: customData,
              name: cat.name,
              meta: cat.name,
              offsetgroup: this.offsetGroups ? this.offsetGroups[sampleName] : sampleName,
              legendgroup: cat.name,
              showlegend: sampleIdx === 0,
              ...traceParams,
              marker: { ...traceParams.marker, color: "rgba(" + cat.color + "," + alpha + ")" },
              hovertemplate: hoverTemplates[cat.name],
            });
          }
        });
      });

      return traces;
    }

    return cats.map((cat) => {
      if (this.layout.barmode !== "group") {
        // Plotting each sample as a separate trace to be able to set alpha for each
        // sample color separately, so we can dim the de-highlighted samples.
        return this.filteredSettings.map((sample, sampleIdx) => {
          let params = JSON.parse(JSON.stringify(traceParams)); // deep copy

          let alpha = highlighted.length > 0 && sample.highlight === null ? 0.1 : 1;
          params.marker.color = "rgba(" + cat.color + "," + alpha + ")";

          return {
            type: "bar",
            x: [cat.data[sampleIdx]],
            y: [sample.name],
            name: cat.name,
            meta: cat.name,
            // To make sure the legend uses bright category colors and not the dim ones:
            showlegend: sampleIdx === firstHighlightedSample,
            legendgroup: cat.name,
            ...params,
          };
        });
      } else {
        // "group"
        // Plotly adds giant gaps between bars in the group mode when adding each sample as a
        // separate trace. Sacrificing dimming the de-highlighted bars to get rid of this gap.
        let params = JSON.parse(JSON.stringify(traceParams)); // deep copy
        let samples = this.filteredSettings.map((s) => s.name);
        params.marker.color = "rgb(" + cat.color + ")";

        return {
          type: "bar",
          x: cat.data,
          y: samples,
          name: cat.name,
          meta: cat.name,
          ...params,
        };
      }
    });
  }

  exportData(format) {
    let [cats] = this.prepData();

    let delim = format === "tsv" ? "\t" : ",";

    let csv = "Sample" + delim + cats.map((cat) => cat.name).join(delim) + "\n";
    for (let i = 0; i < this.filteredSettings.length; i++) {
      csv += this.filteredSettings[i].name + delim + cats.map((cat) => cat.data[i]).join(delim) + "\n";
    }
    return csv;
  }
}

// Make BarPlot globally available
window.BarPlot = BarPlot;
