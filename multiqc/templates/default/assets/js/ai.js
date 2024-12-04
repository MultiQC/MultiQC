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

function formatReportForAi(onlyGeneralStats = false) {
  let result = "";

  // General statistics section
  result += "\n----------------------\n\n";
  result += "Tools used in the report:\n\n";

  Object.entries(aiReportMetadata.tools).forEach(([modAnchor, mod], modIdx) => {
    result += `${modIdx + 1}. ${mod.name}`;
    if (mod.info) result += `\nDescription: ${mod.info}`;
    if (mod.href && mod.href.length > 0) result += `\nLinks: ${mod.href}`;
    if (mod.comment) result += `\nComment: ${mod.comment}`;
    result += "\n\n";
  });

  result += "\n----------------------\n";
  const generalStatsPlot = mqc_plots["general_stats_table"];
  if (generalStatsPlot) {
    result += `\nMultiQC General Statistics (overview of key QC metrics for each sample, across all tools)`;
    result += `\n${generalStatsPlot.prepDataForLlm()}`;
    result += "\n----------------------\n";
  }

  if (!onlyGeneralStats) {
    Object.entries(aiReportMetadata.sections).forEach(([sectionAnchor, section]) => {
      const mod = aiReportMetadata.tools[section.module_anchor];
      result += `\nTool: ${mod.name}`;
      result += "\n" + formatSectionForAi(sectionAnchor);
      result += "\n\n----------------------";
    });
  }

  return result;
}

function formatModuleHeader(mod) {
  let result = `Tool that produced data: ${mod.name}`;
  if (mod.info) result += `\nTool description: ${mod.info}`;
  if (mod.href && mod.href.length > 0) result += `\nTool URL: ${mod.href}`;
  if (mod.comment) result += `\nTool comment: ${mod.comment}`;
  return result;
}

function formatSectionHeader(section) {
  let result = `Section: ${section.name}`;
  if (section.description) result += `\nSection description: ${section.description}`;
  if (section.comment) result += `\nSection comment: ${section.comment}`;
  if (section.helptext) result += `\nSection help text: ${section.helptext}`;
  return result;
}

function formatModAndSectionMetadata(sectionAnchor, moduleAnchor) {
  if (sectionAnchor === "general_stats_table") return "";
  let result = "";
  if (moduleAnchor) {
    const mod = aiReportMetadata.tools[moduleAnchor];
    result += formatModuleHeader(mod) + "\n\n";
  }
  const section = aiReportMetadata.sections[sectionAnchor];
  result += formatSectionHeader(section);
  return result;
}

function formatSectionForAi(sectionAnchor, moduleAnchor, options) {
  if (sectionAnchor === "general_stats_table") {
    return formatReportForAi(true);
  }

  let result = formatModAndSectionMetadata(sectionAnchor, moduleAnchor);
  if (result) result += "\n";

  const section = aiReportMetadata.sections[sectionAnchor];
  if (section.content_before_plot) result += section.content_before_plot + "\n\n";
  if (section.content) result += section.content + "\n\n";

  const plotAnchor = section.plot_anchor;
  let plot = mqc_plots[plotAnchor];
  if (!plot) return result;

  if (plot.pconfig && plot.pconfig.title) result += `Title: ${plot.pconfig.title}\n`;

  result += "\n" + plot.prepDataForLlm(options);

  return result;
}

async function summarizeWithAi(button, options) {
  const { wholeReport, table } = options;

  const sectionAnchor = button.data("section-anchor");
  const moduleAnchor = button.data("module-anchor");

  let content;
  let systemPrompt;
  if (wholeReport) {
    content = formatReportForAi();
    systemPrompt = systemPromptReport;
  } else if (table) {
    content = formatSectionForAi(sectionAnchor, moduleAnchor, { view: "table" });
    systemPrompt = systemPromptPlot;
  } else {
    content = formatSectionForAi(sectionAnchor, moduleAnchor);
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

// Global (report-level) summary generation
async function generateCallback(e) {
  e.preventDefault();
  const button = $(e.currentTarget);

  if (button.hasClass("ai-generate-more-plot")) {
    summarizeWithAi(button, { wholeReport: false, table: false });
  } else if (button.hasClass("ai-generate-more-table")) {
    summarizeWithAi(button, { wholeReport: false, table: true });
  } else {
    summarizeWithAi(button, { wholeReport: false, table: false });
  }
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

  function copyForAi(e, options) {
    e.preventDefault();
    const button = $(e.currentTarget);

    const { wholeReport, table } = options;

    const sectionAnchor = button.data("section-anchor");
    const moduleAnchor = button.data("module-anchor");

    let content;
    let systemPrompt;
    if (wholeReport) {
      content = formatReportForAi();
      systemPrompt = "You are given data of a MultiQC report";
    } else if (table) {
      content = formatSectionForAi(sectionAnchor, moduleAnchor, { view: "table" });
      systemPrompt = "You are given a single MultiQC report table";
    } else {
      content = formatSectionForAi(sectionAnchor, moduleAnchor);
      systemPrompt = "You are given data of a single MultiQC report section with a plot";
    }

    systemPrompt += ". Your task is to analyse the data and give a concise summary.";

    const text = multiqcDescription + "\n" + systemPrompt + "\n\n" + content;

    navigator.clipboard.writeText(text);
    const originalButtonText = button.find(".button-text").text();
    button.find(".button-text").text("Copied!");
    setTimeout(() => {
      button.find(".button-text").text(originalButtonText);
    }, 2000);
  }

  $("button#ai_copy_content_report").click((e) => copyForAi(e, { wholeReport: true, table: false }));
  $("button.ai-copy-content-section").click((e) => copyForAi(e, { wholeReport: false, table: false }));
  $("button.ai-copy-content-table").click((e) => copyForAi(e, { wholeReport: false, table: true }));
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
