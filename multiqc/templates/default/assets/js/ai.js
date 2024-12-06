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

function formatReportForAi(systemPrompt, onlyGeneralStats = false, generalStatsView = "table") {
  let result = "";
  let currentTokens = 0;
  const provider = getStoredProvider();
  const maxTokens =
    {
      OpenAI: 128000,
      Anthropic: 200000,
      "Seqera AI": 128000,
    }[provider] || 128000;

  // Tools section is highest priority - always include
  result += "\n----------------------\n\n";
  result += "Tools used in the report:\n\n";

  Object.entries(aiReportMetadata.tools).forEach(([modAnchor, mod], modIdx) => {
    const toolContent =
      `${modIdx + 1}. ${mod.name}` +
      (mod.info ? `\nDescription: ${mod.info}` : "") +
      (mod.href && mod.href.length > 0 ? `\nLinks: ${mod.href}` : "") +
      (mod.comment ? `\nComment: ${mod.comment}` : "") +
      "\n\n";
    result += toolContent;
  });

  result += "\n----------------------\n";
  currentTokens = estimateTokenCount(systemPrompt + result);

  // General stats are second priority
  const generalStatsPlot = mqc_plots["general_stats_table"];
  if (generalStatsPlot) {
    const statsContent =
      `\nMultiQC General Statistics (overview of key QC metrics for each sample, across all tools)` +
      `\n${generalStatsPlot.formatForAiPrompt(generalStatsView)}` +
      "\n----------------------\n";

    const statsTokens = estimateTokenCount(statsContent);
    if (currentTokens + statsTokens < maxTokens * 0.95) {
      // Leave 5% buffer
      result += statsContent;
      currentTokens += statsTokens;
    } else {
      console.error(
        `Including general stats alone in AI summary would exceed ${provider}'s token limit (${
          currentTokens + statsTokens
        } > ${maxTokens}). Aborting`,
      );
      return result;
    }
  }

  if (!onlyGeneralStats) {
    // Add sections until we approach the token limit
    Object.entries(aiReportMetadata.sections).forEach(([sectionAnchor, section]) => {
      const mod = aiReportMetadata.tools[section.module_anchor];
      let sectionContent = `\nTool: ${mod.name}\n`;
      sectionContent += formatSectionForAi(sectionAnchor);
      sectionContent += "\n\n----------------------";

      const sectionTokens = estimateTokenCount(sectionContent);
      if (currentTokens + sectionTokens < maxTokens * 0.95) {
        result += sectionContent;
        currentTokens += sectionTokens;
      } else {
        console.warn(
          `Truncating report content to fit within ${provider}'s context window. ` +
            `Current tokens: ${currentTokens}, section would add ${sectionTokens} tokens`,
        );
        return result;
      }
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

function formatSectionForAi(sectionAnchor, moduleAnchor, view) {
  if (sectionAnchor === "general_stats_table") {
    return formatReportForAi("", true, view);
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

  result += "\n" + plot.formatForAiPrompt(view);

  return result;
}

async function summarizeWithAi(button, options) {
  const { wholeReport, table } = options;

  const sectionAnchor = button.data("section-anchor") || "global";
  const moduleAnchor = button.data("module-anchor");
  const view = button.data("view");

  const responseDiv = $("#" + sectionAnchor + "_ai_detailed_summary");
  const errorDiv = $("#" + sectionAnchor + "_ai_summary_error");
  const disclaimerDiv = $("#" + sectionAnchor + "_ai_summary_disclaimer");
  const wrapperDiv = $("#" + sectionAnchor + "_ai_summary");

  let content;
  let systemPrompt;
  if (wholeReport) {
    systemPrompt = systemPromptReport;
    content = formatReportForAi(systemPrompt);
  } else if (table) {
    systemPrompt = systemPromptPlot;
    content = formatSectionForAi(sectionAnchor, moduleAnchor, view);
  } else {
    systemPrompt = systemPromptPlot;
    content = formatSectionForAi(sectionAnchor, moduleAnchor, view);
  }

  // Check total tokens before making the request
  const totalTokens = estimateTokenCount(systemPrompt + content);
  const provider = getStoredProvider();
  const maxTokens =
    {
      OpenAI: 128000,
      Anthropic: 200000,
      "Seqera AI": 128000,
    }[provider] || 128000;

  if (totalTokens > maxTokens * 0.95) {
    errorDiv
      .html(
        `Content exceeds ${provider}'s token limit (${totalTokens} > ${maxTokens}). Try analyzing a smaller section.`,
      )
      .show();
    wrapperDiv.show();
    return;
  }

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
    summarizeWithAi(button, { wholeReport: true, table: false });
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
    const view = button.data("view");

    let content;
    let systemPrompt;
    if (wholeReport) {
      systemPrompt = "You are given data of a MultiQC report";
      content = formatReportForAi(systemPrompt);
    } else if (table) {
      systemPrompt = "You are given a single MultiQC report table";
      content = formatSectionForAi(sectionAnchor, moduleAnchor, view);
    } else {
      systemPrompt = "You are given data of a single MultiQC report section with a plot";
      content = formatSectionForAi(sectionAnchor, moduleAnchor, view);
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

function estimateTokenCount(text) {
  const provider = getStoredProvider();

  try {
    if (provider === "OpenAI") {
      // Use GPT-3 Tokenizer if available
      if (window.GPT3Tokenizer) {
        const tokenizer = new window.GPT3Tokenizer({ type: "gpt3" });
        return tokenizer.encode(text).bpe.length;
      }
    } else if (provider === "Anthropic") {
      // Use Claude's tokenizer if available
      if (window.AnthropicTokenizer) {
        return window.AnthropicTokenizer.countTokens(text);
      }
    }
  } catch (e) {
    console.warn("Error using tokenizer:", e);
  }

  // Fallback to character-based estimation
  // Different models have different ratios, but ~4 chars per token is a reasonable estimate
  return Math.ceil(text.length / 4);
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
