////////////////////////////////////////////////
// AI stuff
////////////////////////////////////////////////

window.continueInSeqeraChatHandler = function (event) {
  let el = $(event.currentTarget);
  let seqeraWebsite = el.data("seqera-website");

  // Either report uuid, or encoded system and chat messages
  let generationId = el.data("generation-id");
  let encodedSystemMessage = el.data("encoded-system-message");
  let encodedChatMessages = el.data("encoded-chat-messages");

  let url = seqeraWebsite + "/ask-ai/";
  if (generationId) url += "?generation-id=" + generationId;

  const chatWindow = window.open(url, "_blank");

  // This is used only for Anthropic and OpenAI providers:
  if (encodedSystemMessage && encodedChatMessages) {
    function sendMessage() {
      chatWindow.postMessage(
        {
          type: "chatInitialMessages",
          content: {
            encodedSystemMessage: encodedSystemMessage,
            encodedChatMessages: encodedChatMessages,
          },
        },
        seqeraWebsite,
      );
    }
    setTimeout(sendMessage, 2000);
  }
};

async function formatReportForAi() {
  let contentPrompt = "";

  // General statistics section
  contentPrompt += "\n----------------------\n\n";
  contentPrompt += "Tools used in the report:\n\n";

  Object.entries(aiReportMetadata.modules).forEach(([modAnchor, mod], modIdx) => {
    contentPrompt += `${modIdx + 1}. ${mod.name}`;
    if (mod.info) contentPrompt += `\nDescription: ${mod.info}`;
    if (mod.href && mod.href.length > 0) contentPrompt += `\nURL: ${mod.href}`;
    if (mod.comment) contentPrompt += `\nComment: ${mod.comment}`;
    contentPrompt += "\n\n";
  });

  contentPrompt += "\n----------------------\n";
  const generalStatsPlot = mqc_plots["general_stats_table"];
  if (generalStatsPlot) {
    contentPrompt += `\nMultiQC General Statistics (overview of key QC metrics for each sample, across all tools)`;
    contentPrompt += `\n${generalStatsPlot.prepDataForLlm()}`;
    contentPrompt += "\n----------------------\n";
  }
  // TODO: suffixes and metric descriptions to general stats (and other tables)

  Object.entries(aiReportMetadata.sections).forEach(([sectionAnchor, section]) => {
    let plotContent = section.plot_content;
    let plot = null;
    if (!plotContent) {
      plot = mqc_plots[section.plot_anchor];
      if (plot) {
        if (plot.title) plotContent += `\nTitle: ${plot.title}`;
        plotContent += plot.prepDataForLlm();
      }
    }
    if (!plotContent) return;

    contentPrompt += `\n\nSection: ${section.name}`;
    if (section.description) contentPrompt += `\nDescription: ${section.description}`;
    if (section.comment) contentPrompt += `\nComment: ${section.comment}`;
    if (section.helptext) contentPrompt += `\nHelp text: ${section.helptext}`;

    if (plot) {
      let plotType = plot.plotType;
      if (plotType == "xy_line") plotType = "line";
      if (plotType == "violin") plotType = "table";
      contentPrompt += `\nPlot type: ${plotType}`;
      if (plot.plotDescription) contentPrompt += `\nPlot description: ${plot.plotDescription}`;
      if (plot.plotHelptext) contentPrompt += `\nPlot helptext: ${plot.plotHelptext}`;
    }

    contentPrompt += `\n${plotContent}`;
    contentPrompt += "\n\n----------------------";
  });
  return contentPrompt;
}

function formatSectionMetadata(mod, section) {
  let contentPrompt = "";
  if (mod) {
    let moduleName = mod["name"];
    let moduleSummary = mod["info"];
    let moduleHref = mod["href"];
    let moduleComment = mod["comment"];
    contentPrompt += `\nTool that generated data: ${moduleName}`;
    if (moduleSummary) contentPrompt += `\nTool description: ${moduleSummary}`;
    if (moduleHref) contentPrompt += `\nTool URL: ${moduleHref}`;
    if (moduleComment) contentPrompt += `\nTool comment: ${moduleComment}`;
  }
  if (section) {
    let sectionName = section["name"];
    let sectionDescription = section["description"];
    let sectionComment = section["comment"];
    let sectionHelptext = section["helptext"];
    contentPrompt += `\n\nSection: ${sectionName}`;
    if (sectionDescription) contentPrompt += `\nSection description: ${sectionDescription}`;
    if (sectionComment) contentPrompt += `\nSection comment: ${sectionComment}`;
    if (sectionHelptext) contentPrompt += `\nSection help text: ${sectionHelptext}`;
  }
  return contentPrompt;
}

async function formatSectionForAi(sectionAnchor, moduleAnchor) {
  let result = "";
  let plotAnchor = null;
  if (sectionAnchor !== "general_stats_table") {
    mod = aiReportMetadata.tools[moduleAnchor];
    section = aiReportMetadata.sections[sectionAnchor];
    result += formatSectionMetadata(mod, section);
    plotAnchor = section["plot_anchor"];
  } else {
    plotAnchor = "general_stats_table";
  }

  let plot = mqc_plots[plotAnchor];
  let plotType = plot.plotType;
  if (plotType == "xy_line") plotType = "line";
  if (plotType == "violin") plotType = "table";

  result += `\nPlot type: ${plotType}\n`;
  if (plot.pconfig && plot.pconfig.title) result += `Title: ${plot.pconfig.title}\n`;
  result += "\n" + plot.prepDataForLlm();
  return result;
}

// Global (report-level) summary generation
async function generateCallback(e) {
  e.preventDefault();

  const button = $(e.currentTarget);
  let content;
  let systemPrompt;

  const sectionAnchor = button.data("section-anchor") || "global";
  if (button.data("report-metadata-base64")) {
    const metadata = JSON.parse(atob(button.data("report-metadata-base64")));
    content = await formatReportForAi(metadata);
    systemPrompt = systemPromptReport;
  } else if (button.data("section-anchor") && button.data("module-anchor")) {
    const sectionAnchor = button.data("section-anchor");
    const moduleAnchor = button.data("module-anchor");
    content = await formatSectionForAi(sectionAnchor, moduleAnchor);
    systemPrompt = systemPromptPlot;
  }

  const responseDiv = $("#" + sectionAnchor + "_ai_detailed_summary");
  const errorDiv = $("#" + sectionAnchor + "_ai_summary_error");
  const disclaimerDiv = $("#" + sectionAnchor + "_ai_summary_disclaimer");
  const wrapperDiv = $("#" + sectionAnchor + "_ai_summary");

  let provider = getStoredProvider();
  // Check for stored API key
  let aiApiKey = getStoredApiKey(provider);
  if (!aiApiKey || aiApiKey === undefined) {
    // Open the AI toolbox section
    mqc_toolbox_openclose("#mqc_ai", true);
    return;
  }

  // Disable button and show loading state
  button.prop("disabled", true);
  originalButtonHtml = button.html();
  button.html("Generating...");

  const startTime = performance.now();
  let fullModelName = null;
  await (async () => {
    let receievedMarkdown = "";
    runStreamGeneration({
      systemPrompt: systemPrompt,
      userMessage: content,
      tags: ["multiqc"],
      onStreamStart: (model) => {
        fullModelName = model;
        wrapperDiv.show();
      },
      onStreamNewToken: (token) => {
        receievedMarkdown += token;
        responseDiv.html(markdownToHtml(receievedMarkdown));
      },
      onStreamError: (error) => {
        errorDiv.html(error).show();
        wrapperDiv.show();
        wrapUpResponse(
          button,
          originalButtonHtml,
          responseDiv,
          disclaimerDiv,
          errorDiv,
          wrapperDiv,
          provider,
          getStoredModelName(provider),
        );
      },
      onStreamComplete: () => {
        const provider = getStoredProvider();
        wrapUpResponse(
          button,
          originalButtonHtml,
          responseDiv,
          disclaimerDiv,
          errorDiv,
          wrapperDiv,
          provider,
          fullModelName,
        );
        // Save response to localStorage
        const elementId = button.data("plot-anchor") || "global";
        localStorage.setItem(
          `ai_response_${reportUuid}_${elementId}`,
          JSON.stringify({
            text: receievedMarkdown,
            provider: provider,
            model: fullModelName,
            timestamp: Date.now(),
          }),
        );
        const endTime = performance.now();
        console.log(`Time to generate more: ${endTime - startTime}ms`);
      },
    });
  })();
}

$(function () {
  // "Show More" button to expand pre-generated full AI summary
  $(".ai-summary").each(function () {
    const $details = $(this).find("details");
    const $showMoreBtn = $(this).find(".ai-summary-expand");
    $showMoreBtn.on("click", function (e) {
      if ($details.prop("open")) {
        $details.prop("open", false);
        $showMoreBtn.addClass("ai-summary-expand-closed");
      } else {
        $details.prop("open", true);
        $showMoreBtn.removeClass("ai-summary-expand-closed");
      }
    });

    // Do no expand when clicked on the whole area
    $details.on("click", function (e) {
      e.preventDefault();
    });
  });

  // Click handler for "AI Summary" button to dynamically generate plot summaries
  $("button.ai-generate-more").each(function () {
    const button = $(this);
    const sectionAnchor = button.data("section-anchor") || "global";
    const plotAnchor = button.data("plot-anchor") || "global";

    const responseDiv = $("#" + sectionAnchor + "_ai_detailed_summary").show();
    const disclaimerDiv = $("#" + sectionAnchor + "_ai_summary_disclaimer").show();
    const errorDiv = $("#" + sectionAnchor + "_ai_summary_error");
    const wrapperDiv = $("#" + sectionAnchor + "_ai_summary");

    const originalButtonHtml = button.html();
    const cachedSummaryDump = localStorage.getItem(`ai_response_${reportUuid}_${plotAnchor}`);
    if (cachedSummaryDump) {
      // Load cached AI responses on page load
      const cachedSummary = JSON.parse(cachedSummaryDump);
      responseDiv.html(markdownToHtml(cachedSummary.text));
      if (wrapperDiv) wrapperDiv.show();
      wrapUpResponse(
        button,
        originalButtonHtml,
        responseDiv,
        disclaimerDiv,
        errorDiv,
        wrapperDiv,
        cachedSummary.provider,
        cachedSummary.model,
      );
    } else {
      button.click(generateCallback);
    }
  });

  // Click handler to highlight samples
  $(document).on("click", "sample", function (e) {
    e.preventDefault();
    let sampleName = $(this).text();

    if (sampleName.includes("*")) {
      // Replace * character with regex
      sampleName = sampleName.replace(/\*/g, ".*");
      // Turn on regex mode
      $(".mqc_regex_mode").find(".re_mode").removeClass("off").addClass("on").text("on");
    }

    let color = $(this).css("color");
    let highlightedSamples = window.mqc_highlight_f_texts;
    if (!highlightedSamples.includes(sampleName)) {
      $("#mqc_colour_filter").val(sampleName);
      $("#mqc_colour_filter_color").val(rgbToHex(color));
      $(this).css("font-weight", "bold");
      // also highlight all <sample> elements in text that match the sample name
      $("sample").each(function () {
        if ($(this).text().indexOf(sampleName) > -1) $(this).css("font-weight", "bold");
      });
    } else {
      $("#mqc_col_filters li").each(function () {
        if ($(this).children("input").attr("value") === sampleName) {
          $(this).children(".close").click();
        }
      });
      $(this).css("font-weight", "normal");
      // also remove the bold from all <sample> elements in text that match the sample name
      $("sample").each(function () {
        if ($(this).text().indexOf(sampleName) > -1) $(this).css("font-weight", "normal");
      });
    }

    $("#mqc_color_form").trigger("submit");
    $("#mqc_cols_apply").click();
  });

  // Click handler for "Copy Plot Data" button to copy plot data to clipboard
  $("button.ai-copy-content-plot").click(async function (e) {
    e.preventDefault();

    const button = $(this);
    const sectionAnchor = button.data("section-anchor");
    const moduleAnchor = button.data("module-anchor");
    const plotContent = await formatSectionForAi(sectionAnchor, moduleAnchor);

    const systemPrompt =
      multiqcDescription +
      `\
You are given a single MultiQC report section with a plot. 
Your task is to analyse the data and give a concise summary.
    `;
    navigator.clipboard.writeText(systemPrompt + plotContent);
    const originalButtonText = button.find(".button-text").text();
    button.find(".button-text").text("Copied!");
    setTimeout(() => {
      button.find(".button-text").text(originalButtonText);
    }, 2000);
  });

  // Click handler for "Copy Report Data" button to copy report data to clipboard
  $("button#ai_copy_content_report").click(async function (e) {
    e.preventDefault();

    const systemPrompt =
      multiqcDescription +
      `\
You are given data from such a report. Your task is to analyse this data and generate
a concise summary of the report.
    `;

    const button = $(this);
    const content = await formatReportForAi();

    navigator.clipboard.writeText(systemPrompt + content);
    const originalButtonText = button.find(".button-text").text();
    button.find(".button-text").text("Copied!");
    setTimeout(() => {
      button.find(".button-text").text(originalButtonText);
    }, 2000);
  });

  // Click handler for "Copy Table Data" button to copy table data to clipboard
  $("button.ai-copy-content-table").click(async function (e) {
    e.preventDefault();

    const button = $(this);
    const tableAnchor = button.data("table-anchor");
    const sectionAnchor = button.data("section-anchor");
    const moduleAnchor = button.data("module-anchor");
    let mod = null;
    let section = null;
    if (sectionAnchor != "general_stats_table") {
      mod = aiReportMetadata.tools[moduleAnchor];
      section = aiReportMetadata.sections[sectionAnchor];
    }

    const table = $("#" + tableAnchor);

    // Get table headers
    const headers = [];
    table.find("thead th").each(function () {
      const tt = $(this).find(".mqc_table_tooltip");
      let header;
      if (tt.length > 0) {
        header = tt.data("original-title");
      } else {
        header = $(this).text().trim();
      }
      headers.push(header);
    });

    // Get table data
    const rows = [];
    table.find("tbody tr").each(function () {
      const row = [];
      $(this)
        .find("td,th")
        .each(function () {
          row.push($(this).text().trim());
        });
      if (row.length > 0) rows.push(row);
    });

    // Format data as markdown table
    const content = `
| ${headers.join(" | ")} |
|${headers.map(() => " --- ").join("|")}|
${rows.map((row) => `| ${row.join(" | ")} |`).join("\n")}`;

    const contentPrompt = formatSectionMetadata(mod, section) + content;

    const systemPrompt =
      multiqcDescription +
      `\
You are given a General Statistics table from such a report. Your task is to analyse the data and generate
a concise summary.
    `;
    navigator.clipboard.writeText(systemPrompt + contentPrompt);
    const originalButtonText = button.find(".button-text").text();
    button.find(".button-text").text("Copied!");
    setTimeout(() => {
      button.find(".button-text").text(originalButtonText);
    }, 2000);
  });
});

async function wrapUpResponse(
  button,
  originalButtonHtml,
  responseDiv,
  disclaimerDiv,
  errorDiv,
  wrapperDiv,
  provider,
  model,
) {
  disclaimerDiv.html(`This summary is AI-generated. Provider: ${provider}, model: ${model}`).show();
  const elementId = button.data("plot-anchor") || "global";
  // Change button to "Reset" state
  button
    .text("Clear AI summary")
    .prop("style", "background-color: #f2f2f2;")
    .prop("disabled", false)
    .off("click")
    .on("click", function (e) {
      // Reset and change button back to "Generate" state
      e.preventDefault();
      let sectionAnchor = button.data("section-anchor");
      if (sectionAnchor) {
        $("#" + sectionAnchor + "_ai_summary").hide();
      }
      localStorage.removeItem(`ai_response_${reportUuid}_${elementId}`);
      responseDiv.html("");
      errorDiv.html("");
      disclaimerDiv.html("");
      wrapperDiv.hide();
      button.html(originalButtonHtml).prop("style", "background-color: white;").off("click").click(generateCallback);
    });
}
