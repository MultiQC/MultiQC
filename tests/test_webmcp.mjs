// Run with: node --test tests/test_webmcp.mjs
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { runInNewContext } from "node:vm";

const source = readFileSync(new URL("../multiqc/templates/default/src/js/webmcp.js", import.meta.url), "utf8");

test("WebMCP reads existing plot data, bounds output and only navigates report sections", async () => {
  const registered = new Map();
  const warnings = [];
  const section = {
    scrollIntoView(options) {
      this.scrolled = options;
    },
    setAttribute(name, value) {
      this[name] = value;
    },
    focus(options) {
      this.focused = options;
    },
  };
  const elements = new Map(["general_stats", "qc", "qc_section", "qc_plot", "ai-api-key"].map((id) => [id, section]));
  // Exercise the real table formatter: small scaled counts must retain precision and units.
  const plotWindow = {};
  runInNewContext(
    readFileSync(new URL("../multiqc/templates/default/src/js/plots/violin.js", import.meta.url), "utf8"),
    {
      Plot: class {},
      window: plotWindow,
    },
  );
  const table = Object.create(plotWindow.ViolinPlot.prototype);
  table.pconfig = { col1_header: "Sample" };
  table.prepData = () => [
    ["reads"],
    { reads: { title: "Reads", description: "Read count", suffix: " M" } },
    ["Sample A", "Sample B"],
    ["Sample A", "Sample B"].map((name) => ({ name, originalName: name })),
    { reads: { "Sample A": 0.0001, "Sample B": 1 } },
    { reads: {} },
  ];
  let text = table.formatDatasetForAiPrompt({});
  assert.match(text, /\|Sample A\|0\.0001 M\|/);
  assert.match(text, /\|Sample B\|1 M\|/);
  const plot = {
    anchor: "qc_plot",
    pconfig: { title: "QC metrics" },
    plotType: "violin plot",
    isDownsampled: true,
    formatForAiPrompt: () => text,
  };
  const window = {
    callAfterDecompressed: [],
    aiPseudonymMap: { "Sample A": "SAMPLE_1", "Sample B": "SAMPLE_2" },
    aiReportMetadata: {
      tools: { qc: { name: "QC tool" }, hidden: { name: "Hidden module" } },
      sections: {
        qc_section: { name: "QC section", module_anchor: "qc", plot_anchor: "qc_plot" },
        hidden_section: { name: "Hidden section" },
      },
    },
    mqc_plots: { qc_plot: plot, unknown_type: null, absent: { anchor: "not_in_report" } },
  };
  const document = {
    title: "Test MultiQC report",
    getElementById: (id) => elements.get(id),
  };
  const console = { warn: (...args) => warnings.push(args) };
  runInNewContext(source, { window, document, console });
  assert.equal(window.callAfterDecompressed.length, 1);
  const initialize = window.callAfterDecompressed[0];
  await initialize(); // Unsupported browsers leave the report alone.
  assert.equal(registered.size, 0);
  document.modelContext = {
    registerTool: async (tool) => registered.set(tool.name, tool),
  };
  await initialize();
  assert.equal(registered.size, 3);
  for (const tool of registered.values()) {
    assert.equal(tool.inputSchema.type, "object");
    assert.equal(tool.annotations.untrustedContentHint, true);
  }
  const summary = await registered.get("multiqc_get_report_summary").execute({});
  assert.equal(summary.sample_count, 2);
  assert.equal(summary.modules.length, 1);
  assert.equal(summary.sections.length, 2);
  assert.equal(summary.plots.length, 1);
  assert.equal(summary.plots[0].anchor, "qc_plot");
  const read = registered.get("multiqc_get_plot_data").execute;
  assert.equal((await read({ anchor: "qc_plot" })).text, text);
  assert.equal((await read({ anchor: "qc_plot" })).downsampled, true);
  let reconstructed = "";
  let offset = 0;
  do {
    const page = await read({ anchor: "qc_plot", offset, limit: 5 });
    reconstructed += page.text;
    offset = page.next_offset;
  } while (offset !== null);
  assert.equal(reconstructed, text);
  assert.equal((await read({ anchor: "qc_plot", offset: 500 })).text, "");
  for (const anchor of ["missing", "toString", "__proto__", "unknown_type", "absent", null]) {
    assert.ok((await read({ anchor })).error);
  }
  for (const args of [{ offset: -1 }, { limit: 0 }, { limit: 50001 }, { limit: 1.5 }, { offset: "0" }]) {
    assert.ok((await read({ anchor: "qc_plot", ...args })).error);
  }
  text = "SAMPLE_1: QC 12.34%"; // The existing formatter controls filters and anonymization.
  assert.equal((await read({ anchor: "qc_plot" })).text, text);
  const show = registered.get("multiqc_show_section");
  assert.equal(show.annotations.readOnlyHint, false);
  assert.ok((await show.execute({ anchor: "ai-api-key" })).error);
  assert.equal(section.scrolled, undefined);
  assert.equal((await show.execute({ anchor: "qc_section" })).anchor, "qc_section");
  assert.equal(section.scrolled.block, "start");
  assert.equal(section.tabindex, "-1");
  assert.equal(section.focused.preventScroll, true);
  document.modelContext.registerTool = async () => {
    throw new Error("Blocked by permissions policy");
  };
  await initialize();
  assert.equal(warnings.length, 3); // Registration failures never reject report initialization.
});
