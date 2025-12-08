class BarPlot extends Plot {
  constructor(dump) {
    super(dump);
    this.filteredSettings = [];
    this.groupSettingsMap = null; // For multicategory: maps original group name -> toolbox settings
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

    // Process group labels through toolbox if present (for hiding/renaming/highlighting groups)
    if (groupLabels && groupLabels.length > 0) {
      let uniqueGroups = [...new Set(groupLabels)];
      let groupSettings = applyToolboxSettings(uniqueGroups);
      this.groupSettingsMap = Object.fromEntries(uniqueGroups.map((g, i) => [g, groupSettings[i]]));

      // Determine which rows to keep based on group visibility only (not sample visibility)
      let keepRow = groupLabels.map((gl) => {
        let groupSetting = this.groupSettingsMap[gl];
        return !(groupSetting && groupSetting.hidden);
      });

      // For multicategory, filteredSettings just tracks the samples (no toolbox ops on samples)
      this.filteredSettings = samples.filter((_, si) => keepRow[si]).map((name) => ({ name }));

      // Store original group labels (for data indexing) and renamed labels (for display)
      this.originalGroupLabels = groupLabels.filter((_, si) => keepRow[si]);
      this.filteredGroupLabels = this.originalGroupLabels.map((gl) => this.groupSettingsMap[gl]?.name || gl);

      this.offsetGroups = offsetGroups;

      cats = cats.map((cat) => {
        let data = this.pActive ? cat["data_pct"] : cat.data;
        return {
          data: data.filter((_, si) => keepRow[si]),
          color: cat.color,
          name: cat.name,
        };
      });

      return [cats];
    }

    // Non-multicategory: existing logic with full sample toolbox support
    this.groupSettingsMap = null;
    this.originalGroupLabels = null;

    // Rename and filter samples:
    this.filteredSettings = samplesSettings.filter((s) => !s.hidden);

    // Filter group labels to match visible samples
    this.filteredGroupLabels = null;
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
      // Check if any group highlights are active
      let anyGroupHighlights = Object.values(this.groupSettingsMap || {}).some((g) => g.highlight);

      // Build data structures using ORIGINAL group labels as keys (not renamed)
      let uniqueSamples = [];
      let seenSamples = new Set();
      let sampleGroupData = {}; // sampleName -> originalGroupLabel -> [{name, value}, ...]
      let sampleGroupEntries = {}; // sampleName -> [{originalGroupLabel, displayGroupLabel, dataIdx}, ...]

      this.filteredSettings.forEach((sample, idx) => {
        let sampleName = sample.name;
        let originalGroupLabel = this.originalGroupLabels[idx];
        let displayGroupLabel = this.filteredGroupLabels[idx];

        // Track unique samples
        if (!seenSamples.has(sampleName)) {
          seenSamples.add(sampleName);
          uniqueSamples.push(sampleName);
          sampleGroupData[sampleName] = {};
          sampleGroupEntries[sampleName] = [];
        }

        // Store group entry for this sample (use original label as key)
        sampleGroupEntries[sampleName].push({ originalGroupLabel, displayGroupLabel, dataIdx: idx });

        // Build category data for this sample+group combination (use original label as key)
        if (!sampleGroupData[sampleName][originalGroupLabel]) {
          sampleGroupData[sampleName][originalGroupLabel] = cats.map((cat) => ({
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
        cats.forEach((cat) => {
          // Collect data for this sample across all its groups
          let displayGroupLabels = [];
          let groupData = [];
          let customData = [];

          sampleGroupEntries[sampleName].forEach(({ originalGroupLabel, displayGroupLabel, dataIdx }) => {
            // Check if this group is highlighted
            let groupHighlight = this.groupSettingsMap[originalGroupLabel]?.highlight;
            let alpha = anyGroupHighlights && !groupHighlight ? 0.1 : 1;

            displayGroupLabels.push(displayGroupLabel);
            groupData.push(cat.data[dataIdx]);
            customData.push({
              sampleName: sampleName,
              catData: sampleGroupData[sampleName][originalGroupLabel],
              alpha: alpha,
            });
          });

          if (displayGroupLabels.length > 0) {
            // Calculate alpha per bar based on group highlight
            let alphas = customData.map((d) => d.alpha);
            let colors = alphas.map((a) => "rgba(" + cat.color + "," + a + ")");

            traces.push({
              type: "bar",
              x: groupData,
              y: displayGroupLabels,
              customdata: customData,
              name: cat.name,
              meta: cat.name,
              offsetgroup: this.offsetGroups ? this.offsetGroups[sampleName] : sampleName,
              legendgroup: cat.name,
              showlegend: sampleIdx === 0,
              ...traceParams,
              marker: { ...traceParams.marker, color: colors },
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

  afterPlotCreated() {
    super.afterPlotCreated();

    // Only apply in multicategory mode with group highlights
    if (!this.groupSettingsMap) return;

    const plotDiv = document.getElementById(this.anchor);
    if (!plotDiv) return;

    // Style Y-axis tick labels for highlighted groups
    plotDiv.querySelectorAll(".yaxislayer-above .ytick text, .yaxislayer-above .ytick tspan").forEach((tick) => {
      const groupName = tick.textContent;
      // Find by renamed name, get original's highlight
      let origGroup = Object.keys(this.groupSettingsMap).find((k) => this.groupSettingsMap[k].name === groupName);
      let setting = this.groupSettingsMap[origGroup];
      if (setting?.highlight) {
        tick.style.fill = setting.highlight;
        tick.style.fontWeight = "bold";
      }
    });
  }

  exportData(format) {
    let [cats] = this.prepData();

    let delim = format === "tsv" ? "\t" : ",";

    // For multicategory, include group column (use display labels which may be renamed)
    let useMulticategory = this.filteredGroupLabels && this.filteredGroupLabels.length > 0;

    if (useMulticategory) {
      let csv = "Group" + delim + "Sample" + delim + cats.map((c) => c.name).join(delim) + "\n";
      for (let i = 0; i < this.filteredSettings.length; i++) {
        csv +=
          this.filteredGroupLabels[i] +
          delim +
          this.filteredSettings[i].name +
          delim +
          cats.map((c) => c.data[i]).join(delim) +
          "\n";
      }
      return csv;
    }

    // Non-multicategory export
    let csv = "Sample" + delim + cats.map((cat) => cat.name).join(delim) + "\n";
    for (let i = 0; i < this.filteredSettings.length; i++) {
      csv += this.filteredSettings[i].name + delim + cats.map((cat) => cat.data[i]).join(delim) + "\n";
    }
    return csv;
  }
}

// Make BarPlot globally available
window.BarPlot = BarPlot;
