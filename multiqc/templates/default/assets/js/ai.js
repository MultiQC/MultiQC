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

function formatReportForAi(systemTokens, onlyGeneralStats = false, generalStatsView = "table") {
  let userPrompt = "";
  let currentTokens = systemTokens ?? 0;
  const providerId = $("#ai-provider").val();
  const provider = AI_PROVIDERS[providerId];
  const model = $("#ai-model").val();

  const maxTokens = getMaxTokens(model);

  // Tools section is highest priority - always include
  userPrompt += "\n----------------------\n\n";
  userPrompt += "Tools used in the report:\n\n";

  Object.entries(aiReportMetadata.tools).forEach(([modAnchor, mod], modIdx) => {
    const toolsContext =
      `${modIdx + 1}. ${mod.name}` +
      (mod.info ? `\nDescription: ${mod.info}` : "") +
      (mod.href && mod.href.length > 0 ? `\nLinks: ${mod.href}` : "") +
      (mod.comment ? `\nComment: ${mod.comment}` : "") +
      "\n\n";
    userPrompt += toolsContext;
  });

  userPrompt += "\n----------------------\n";
  currentTokens += countTokens(userPrompt);

  // General stats are second priority
  const generalStatsPlot = mqc_plots["general_stats_table"];
  if (generalStatsPlot) {
    const genStatsContext =
      `\nMultiQC General Statistics (overview of key QC metrics for each sample, across all tools)` +
      `\n${generalStatsPlot.formatForAiPrompt(generalStatsView)}` +
      "\n----------------------\n";

    userPrompt += genStatsContext;
    const genStatsTokens = countTokens(genStatsContext);
    if (currentTokens + genStatsTokens <= maxTokens) {
      currentTokens += genStatsTokens;
    } else {
      console.error(
        `General stats alone would already exceed the token limit of ${provider.name} (${
          currentTokens + genStatsTokens
        } > ${maxTokens}). Cannot summarize the report`,
      );
      return userPrompt;
    }
  }

  let sectionsContext = "";
  if (!onlyGeneralStats) {
    // Add sections until we approach the token limit
    for (const [sectionAnchor, section] of Object.entries(aiReportMetadata.sections)) {
      const mod = aiReportMetadata.tools[section.module_anchor];
      sectionsContext += `\nTool: ${mod.name}\n`;
      sectionsContext += formatSectionForAi(sectionAnchor);
      sectionsContext += "\n\n----------------------";

      const sectionsTokens = countTokens(sectionsContext);
      if (currentTokens + sectionsTokens > maxTokens) {
        console.warn(
          `Truncating prompt to only the general stats to fit within the context window of ${provider.name} (${maxTokens} tokens). ` +
            `Tokens estimate: ${currentTokens}, with sections: at least ${currentTokens + sectionsTokens}`,
        );
        return userPrompt; // Stop iterating through sections and return only general stats
      }
    }
  }
  userPrompt += sectionsContext;
  return userPrompt;
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

function formatSectionForAi(sectionAnchor, moduleAnchor, plotView) {
  let result = formatModAndSectionMetadata(sectionAnchor, moduleAnchor);
  if (result) result += "\n";

  const section = aiReportMetadata.sections[sectionAnchor];
  if (section.content_before_plot) result += section.content_before_plot + "\n\n";
  if (section.content) result += section.content + "\n\n";

  const plotAnchor = section.plot_anchor;
  let plot = mqc_plots[plotAnchor];
  if (!plot) return result;

  if (plot.pconfig && plot.pconfig.title) result += `Title: ${plot.pconfig.title}\n`;

  result += "\n" + plot.formatForAiPrompt(plotView);

  return result;
}

function getMaxTokens(model) {
  if (model.startsWith("gpt")) return 128000;
  if (model.startsWith("claude")) return 200000;
  return 128000;
}

async function summarizeWithAi(button) {
  const isGlobal = button.hasClass("ai-generate-button-global");
  const isMore = button.hasClass("ai-generate-button-more");

  const responseDiv = $("#" + button.data("response-div"));
  const errorDiv = $("#" + button.data("error-div")).hide();
  const disclaimerDiv = $("#" + button.data("disclaimer-div"));
  const wrapperDiv = $("#" + button.data("wrapper-div"));

  const sectionAnchor = button.data("section-anchor") || "global";
  const moduleAnchor = button.data("module-anchor");
  const plotView = button.data("plot-view");
  const clearText = button.data("clear-text");

  let content;
  let systemPrompt;
  if (isGlobal) {
    systemPrompt = isMore ? systemPromptReportFull : systemPromptReportShort;
    content = formatReportForAi(countTokens(systemPrompt));
  } else if (sectionAnchor === "general_stats_table") {
    systemPrompt = systemPromptPlot;
    content = formatReportForAi(countTokens(systemPrompt), true, plotView);
  } else {
    systemPrompt = systemPromptPlot;
    content = formatSectionForAi(sectionAnchor, moduleAnchor, plotView);
  }

  // Check total tokens before making the request
  const totalTokens = countTokens(systemPrompt + content);
  const providerId = $("#ai-provider").val();
  const provider = AI_PROVIDERS[providerId];
  let modelName = $("#ai-model").val();
  if (!modelName) {
    modelName = provider.defaultModel;
    $("#ai-model").val(modelName);
    storeModelName(providerId, modelName);
  }
  const maxTokens = getMaxTokens(modelName);

  if (totalTokens > maxTokens) {
    errorDiv.html(`Content exceeds the token limit of ${provider.name} (${totalTokens} > ${maxTokens})`).show();
    if (wrapperDiv) wrapperDiv.show();
    return;
  }

  // Check for stored API key
  let aiApiKey = $("#ai-api-key").val();
  if (!aiApiKey || aiApiKey === undefined) {
    // Open the AI toolbox section
    mqc_toolbox_openclose("#mqc_ai", true);
    return;
  }

  // Disable button and show loading state
  button.prop("disabled", true).html(`Requesting ${provider.name}...`);

  function wrapUpResponse() {
    disclaimerDiv.html(`This summary is AI-generated. Provider: ${provider.name}, model: ${modelName}`).show();
    button.data("action", "clear").prop("disabled", false).html(clearText).addClass("ai-local-content");
  }

  const startTime = performance.now();
  await (async () => {
    let receievedMarkdown = "";
    runStreamGeneration({
      systemPrompt: systemPrompt,
      userMessage: content,
      tags: ["multiqc"],
      onStreamStart: (resolvedModelName) => {
        modelName = resolvedModelName;
        if (wrapperDiv) wrapperDiv.show();
        responseDiv.show();
        button.html(`Generating...`);
      },
      onStreamNewToken: (token) => {
        receievedMarkdown += token;
        responseDiv.html(markdownToHtml(receievedMarkdown));
      },
      onStreamError: (error) => {
        errorDiv.html(error).show();
        wrapUpResponse();
        if (wrapperDiv) wrapperDiv.show();
        if (!isMore && isGlobal) $("#global_ai_summary_more_button_and_disclaimer").hide();
        disclaimerDiv.hide();
      },
      onStreamComplete: () => {
        wrapUpResponse();
        // Save response to localStorage
        const elementId = button.data("plot-anchor") || "global";
        localStorage.setItem(
          `ai_response_${reportUuid}_${elementId}${isMore ? "_more" : ""}`,
          JSON.stringify({
            text: receievedMarkdown,
            provider: providerId,
            model: modelName,
            timestamp: Date.now(),
          }),
        );
        const endTime = performance.now();
        console.log(`Time to generate more: ${endTime - startTime}ms`);
        if (!isMore && isGlobal) $("#global_ai_summary_more_button_and_disclaimer").show();
        disclaimerDiv.show();
      },
    });
  })();
}

// Global (report-level) summary generation
async function generateCallback(e) {
  e.preventDefault();
  const button = $(e.currentTarget);
  const isMore = button.hasClass("ai-generate-button-more");
  const action = button.data("action");
  const responseDiv = $("#" + button.data("response-div"));
  const errorDiv = $("#" + button.data("error-div"));
  const wrapperDiv = $("#" + button.data("wrapper-div"));
  const originalButtonHtml = button.data("original-html");
  const elementId = button.data("plot-anchor") || "global";

  if (action === "clear") {
    e.preventDefault();
    localStorage.removeItem(`ai_response_${reportUuid}_${elementId}${isMore ? "_more" : ""}`);
    responseDiv.html("").hide();
    errorDiv.html("").hide();
    if (wrapperDiv) wrapperDiv.hide();
    button.html(originalButtonHtml).data("action", "generate").removeClass("ai-local-content");
  } else {
    summarizeWithAi(button);
  }
}

$(function () {
  $("#global_ai_summary_expand").each(function () {
    const responseDiv = $("#global_ai_summary_detailed_analysis_response");
    const expandBtn = $("#global_ai_summary_expand");
    const expandBtnGlyphicon = expandBtn.find(".glyphicon");

    let isExpanded = $(this).hasClass("ai-summary-expand-expanded");
    const storedState = localStorage.getItem("mqc_ai_global_summary_expanded");
    if (storedState === "expanded") isExpanded = true;
    if (storedState === "collapsed") isExpanded = false;

    if (isExpanded) {
      responseDiv.show();
      expandBtn.addClass("ai-summary-expand-expanded");
      expandBtnGlyphicon.addClass("glyphicon-chevron-up");
      expandBtnGlyphicon.removeClass("glyphicon-chevron-down");
    } else {
      responseDiv.hide();
      expandBtn.removeClass("ai-summary-expand-expanded");
      expandBtnGlyphicon.addClass("glyphicon-chevron-down");
      expandBtnGlyphicon.removeClass("glyphicon-chevron-up");
    }

    expandBtn.on("click", (e) => {
      e.preventDefault();
      isExpanded = !isExpanded;
      if (isExpanded) {
        responseDiv.show();
        expandBtn.addClass("ai-summary-expand-expanded");
      } else {
        responseDiv.hide();
        expandBtn.removeClass("ai-summary-expand-expanded");
      }
      localStorage.setItem("mqc_ai_global_summary_expanded", isExpanded ? "expanded" : "collapsed");
    });

    // Do no expand when clicked on the whole area
    responseDiv.on("click", function (e) {
      e.preventDefault();
    });
  });

  // Click handler for "AI Summary" button to dynamically generate plot summaries
  $("button.ai-generate-button").each(function () {
    const button = $(this);

    const originalButtonHtml = button.html();
    button.data("original-html", originalButtonHtml).removeClass("ai-local-content");
    const clearText = button.data("clear-text");

    const isMore = button.hasClass("ai-generate-button-more");

    const action = button.data("action");

    const responseDiv = $("#" + button.data("response-div"));
    const disclaimerDiv = $("#" + button.data("disclaimer-div"));
    const wrapperDiv = $("#" + button.data("wrapper-div"));
    if (wrapperDiv) wrapperDiv.addClass("ai-local-content");

    if (action === "clear") {
      button.html(clearText).addClass("ai-local-content");
    } else {
      // Load cached AI responses on page load
      const plotAnchor = button.data("plot-anchor") || "global";
      const cachedSummaryDump = localStorage.getItem(`ai_response_${reportUuid}_${plotAnchor}${isMore ? "_more" : ""}`);
      if (cachedSummaryDump) {
        const cachedSummary = JSON.parse(cachedSummaryDump);
        responseDiv.show().html(markdownToHtml(cachedSummary.text));
        if (wrapperDiv) wrapperDiv.show();
        const provider = AI_PROVIDERS[cachedSummary.provider];
        disclaimerDiv
          .html(`This summary is AI-generated. Provider: ${provider.name}, model: ${cachedSummary.model}`)
          .show();
        button.html(clearText).data("action", "clear").prop("disabled", false).addClass("ai-local-content");
      }
    }

    button.on("click", generateCallback);
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
    const plotView = button.data("plot-view");

    let content;
    let systemPrompt;
    if (wholeReport) {
      systemPrompt = "You are given data of a MultiQC report";
      content = formatReportForAi();
    } else if (sectionAnchor === "general_stats_table") {
      systemPrompt = "You are given the general statistics report table";
      content = formatReportForAi(0, true, plotView);
    } else if (table) {
      systemPrompt = "You are given a single MultiQC report table";
      content = formatSectionForAi(sectionAnchor, moduleAnchor, plotView);
    } else {
      systemPrompt = "You are given data of a single MultiQC report section with a plot";
      content = formatSectionForAi(sectionAnchor, moduleAnchor, plotView);
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

  $("button.ai-copy-content-report").click((e) => copyForAi(e, { wholeReport: true, table: false }));
  $("button.ai-copy-content-plot").click((e) => copyForAi(e, { wholeReport: false, table: false }));
  $("button.ai-copy-content-table").click((e) => copyForAi(e, { wholeReport: false, table: true }));
});

function countTokens(text) {
  // Use simple character-based estimation. Since the data has a lot of numbers
  // and repetitive special characters, using a more conservative estimate than
  // 4 chars per token suitable for real texts - 1.5 chars per token
  return Math.ceil(text.length / 1.5);
}

// Load tokenizers if available
document.addEventListener("DOMContentLoaded", function () {
  // Load GPT-3 Tokenizer
  const gpt3Script = document.createElement("script");
  gpt3Script.src = "https://cdn.jsdelivr.net/npm/gpt3-tokenizer/dist/gpt3-tokenizer.min.js";
  gpt3Script.async = true;

  // Load Claude's Tokenizer
  const claudeScript = document.createElement("script");
  claudeScript.src = "https://cdn.jsdelivr.net/npm/@anthropic-ai/tokenizer/dist/tokenizer.min.js";
  claudeScript.async = true;

  document.head.appendChild(gpt3Script);
  document.head.appendChild(claudeScript);
});
