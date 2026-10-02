// Register after plot initialization, including plots whose rendering is deferred.
window.callAfterDecompressed.push(async function () {
  const statusElements = document.querySelectorAll("[data-webmcp-status]");
  if (typeof document.modelContext?.registerTool !== "function") {
    for (const element of statusElements) element.textContent = "Unavailable in this browser";
    return;
  }

  const sections = Object.entries(window.aiReportMetadata.sections)
    .filter(([anchor]) => document.getElementById(anchor))
    .map(([anchor, section]) => ({
      anchor,
      name: section.name,
      module_anchor: section.module_anchor,
      plot_anchor: section.plot_anchor,
    }));
  if (document.getElementById("general_stats")) {
    sections.unshift({ anchor: "general_stats", name: "General Statistics", plot_anchor: "general_stats_table" });
  }

  const tools = [
    {
      name: "multiqc_get_report_summary",
      description:
        "List this MultiQC report's modules, sections and interactive plot anchors before requesting QC data.",
      inputSchema: { type: "object", properties: {}, additionalProperties: false },
      annotations: { readOnlyHint: true, untrustedContentHint: true },
      execute: async () => ({
        title: document.title,
        sample_count: Object.keys(window.aiPseudonymMap).length,
        modules: Object.entries(window.aiReportMetadata.tools)
          .filter(([anchor]) => document.getElementById(anchor))
          .map(([anchor, module]) => ({ anchor, name: module.name })),
        sections,
        plots: Object.values(window.mqc_plots)
          .filter((plot) => plot && (document.getElementById(plot.anchor) || document.getElementById(plot.tableAnchor)))
          .map((plot) => ({ anchor: plot.anchor, title: plot.pconfig.title, plot_type: plot.plotType })),
      }),
    },
    {
      name: "multiqc_get_plot_data",
      description:
        "Read QC data from an interactive plot, including General Statistics. Uses the current sample/column filters " +
        "and AI sample anonymization. Values can be scaled or downsampled. Returns text with character offsets for paging.",
      inputSchema: {
        type: "object",
        properties: {
          anchor: { type: "string", minLength: 1, description: "A plot anchor from multiqc_get_report_summary." },
          offset: {
            type: "integer",
            minimum: 0,
            default: 0,
            description: "Character offset, or next_offset from a previous result.",
          },
          limit: {
            type: "integer",
            minimum: 1,
            maximum: 50000,
            default: 10000,
            description: "Maximum characters to return in this page.",
          },
        },
        required: ["anchor"],
        additionalProperties: false,
      },
      annotations: { readOnlyHint: true, untrustedContentHint: true },
      execute: async ({ anchor, offset = 0, limit = 10000 }) => {
        if (typeof anchor !== "string" || !Object.hasOwn(window.mqc_plots, anchor) || !window.mqc_plots[anchor]) {
          return { error: "Unknown plot anchor. Use multiqc_get_report_summary to find an interactive plot." };
        }
        if (!Number.isSafeInteger(offset) || offset < 0 || !Number.isSafeInteger(limit) || limit < 1 || limit > 50000) {
          return { error: "offset must be a nonnegative integer and limit must be an integer between 1 and 50000." };
        }
        const plot = window.mqc_plots[anchor];
        if (!document.getElementById(plot.anchor) && !document.getElementById(plot.tableAnchor)) {
          return { error: "This plot is not included in the displayed report." };
        }
        // ponytail: format the whole plot before paging; page datasets first if large-report profiling warrants it.
        const text = plot.formatForAiPrompt(anchor === "general_stats_table" ? "table" : undefined);
        const end = Math.min(offset + limit, text.length);
        return {
          anchor,
          title: plot.pconfig.title,
          plot_type: plot.plotType,
          downsampled: Boolean(plot.isDownsampled),
          text: text.slice(offset, end),
          total_characters: text.length,
          next_offset: end < text.length ? end : null,
        };
      },
    },
    {
      name: "multiqc_show_section",
      description: "Scroll to a section of this MultiQC report. Changes the report view without changing QC data.",
      inputSchema: {
        type: "object",
        properties: {
          anchor: { type: "string", minLength: 1, description: "A section anchor from the report summary." },
        },
        required: ["anchor"],
        additionalProperties: false,
      },
      annotations: { readOnlyHint: false, consequentialHint: false, untrustedContentHint: true },
      execute: async ({ anchor }) => {
        if (!sections.some((section) => section.anchor === anchor)) {
          return { error: "Unknown section anchor. Use multiqc_get_report_summary to find a report section." };
        }
        const target = document.getElementById(anchor);
        target.scrollIntoView({ behavior: "instant", block: "start" });
        // Keep keyboard and assistive technology focus in sync with the visible section.
        target.setAttribute("tabindex", "-1");
        target.focus({ preventScroll: true });
        return { anchor };
      },
    },
  ];

  let registered = 0;
  for (const tool of tools) {
    try {
      await document.modelContext.registerTool(tool);
      registered++;
    } catch (error) {
      // Permission policies or experimental API changes must not break the report.
      console.warn(`Could not register WebMCP tool ${tool.name}:`, error);
    }
  }
  for (const element of statusElements) {
    element.textContent =
      registered === tools.length
        ? `${registered} tools available`
        : registered > 0
          ? `${registered} of ${tools.length} tools available`
          : "Tools could not be enabled";
  }
});
