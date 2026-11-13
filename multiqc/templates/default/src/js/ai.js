////////////////////////////////////////////////
// AI stuff
////////////////////////////////////////////////

import {
  claudeDefaultMaxTokens,
  FALLBACK_CONTEXT_WINDOW,
  parseContextOverflow,
  parseModel,
  ruleContextWindow,
  shrinkBudget,
} from "./ai-models.js";

window.continueInSeqeraChatHandler = function (event) {
  let el = $(event.currentTarget);
  let seqeraWebsite = el.data("seqera-website");

  // Either report uuid, or encoded system and chat messages
  let threadId = el.data("thread-id");

  let url = seqeraWebsite;
  if (threadId) url += "/chat/" + threadId;

  window.open(url, "_blank");
};

// Fit as much of the report as possible into the token budget: tools and General Statistics first,
// then every section that still fits, in report order. Sections that don't fit are listed by name.
function formatReportForAi(systemTokens, onlyGeneralStats = false, generalStatsView = "table", budget = null) {
  let userPrompt = "";
  let currentTokens = systemTokens ?? 0;
  const model = $("#ai-model").val();
  const maxTokens = budget ?? getPromptBudget(model);

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
    currentTokens += countTokens(genStatsContext);
    if (currentTokens > maxTokens) {
      console.error(
        `General stats alone exceed the token budget for ${model} (${currentTokens} > ${maxTokens}). ` +
          "Cannot summarize the report",
      );
      return userPrompt;
    }
  }

  if (!onlyGeneralStats) {
    const skipped = [];
    for (const [sectionAnchor, section] of Object.entries(aiReportMetadata.sections)) {
      const mod = aiReportMetadata.tools[section.module_anchor];
      const sectionContext = `\nTool: ${mod.name}\n` + formatSectionForAi(sectionAnchor) + "\n\n----------------------";
      const sectionTokens = countTokens(sectionContext);
      if (currentTokens + sectionTokens > maxTokens) {
        skipped.push(`${mod.name}: ${section.name}`);
        continue;
      }
      userPrompt += sectionContext;
      currentTokens += sectionTokens;
    }
    if (skipped.length > 0) {
      console.warn(
        `Left out ${skipped.length} report sections to fit the context window of ${model} (budget ${maxTokens} tokens)`,
      );
      userPrompt +=
        "\n----------------------\n\nThese report sections were left out to fit the model's context window: " +
        skipped.join("; ") +
        "\n";
    }
  }
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

// The model's context window: the configured override, then what MultiQC found out when building
// the report, then the model family rules, then a conservative default
function getContextWindow(model) {
  const providerId = $("#ai-provider").val();
  if (providerId === "seqera") return 200000;
  if (providerId === "custom") {
    const custom = parseInt($("#ai-context-window").val());
    if (custom > 0) return custom;
  } else if (aiConfigCustomContextWindow && aiConfigCustomContextWindow !== "None") {
    const configured = parseInt(aiConfigCustomContextWindow);
    if (configured > 0) return configured;
  }
  const config = window.aiModelConfig || {};
  const known = (config.contextWindows || {})[model];
  if (known) return known;
  return ruleContextWindow(parseModel(model)) || config.fallbackContextWindow || FALLBACK_CONTEXT_WINDOW;
}

// Tokens available for the prompt, leaving room for the answer
function getPromptBudget(model) {
  const config = window.aiModelConfig || {};
  const info = parseModel(model);
  const outputTokens =
    config.maxCompletionTokens || (info.vendor === "anthropic" ? claudeDefaultMaxTokens(info) : 4096);
  return Math.max(1, Math.floor(getContextWindow(model) * 0.95) - outputTokens);
}

async function summarizeWithAi(button) {
  const isGlobal = button.hasClass("ai-generate-button-global");
  const isMore = button.hasClass("ai-generate-button-more");

  const responseDiv = $("#" + button.data("response-div"));
  const errorDiv = $("#" + button.data("error-div")).hide();
  const disclaimerDiv = $("#" + button.data("disclaimer-div"));
  const wrapperDiv = $("#" + button.data("wrapper-div"));
  const continueInChatButton = $("#" + button.data("continue-in-chat-button"));

  const sectionAnchor = button.data("section-anchor") || "global";
  const moduleAnchor = button.data("module-anchor");
  const plotView = button.data("plot-view");
  const clearText = button.data("clear-text");

  let title = configTitle == "None" ? "" : configTitle + ": ";

  let systemPrompt;
  // Builds the report data, fitting it into the budget
  let buildContent;
  if (isGlobal) {
    systemPrompt = isMore ? window.systemPromptReportFull : window.systemPromptReportShort;
    buildContent = (budget) => formatReportForAi(countTokens(systemPrompt), false, "table", budget);
    title += "MultiQC report";
  } else if (sectionAnchor === "general_stats_table") {
    systemPrompt = window.systemPromptPlot;
    buildContent = (budget) => formatReportForAi(countTokens(systemPrompt), true, plotView, budget);
    title += "MultiQC General Statistics";
  } else {
    systemPrompt = window.systemPromptPlot;
    buildContent = () => formatSectionForAi(sectionAnchor, moduleAnchor, plotView);
    const section = aiReportMetadata.sections[sectionAnchor];
    title += `MultiQC ${section.name}`;
  }

  const providerId = $("#ai-provider").val();
  const provider = window.AI_PROVIDERS[providerId];

  let modelName = $("#ai-model").val();
  let aiApiKey = $("#ai-api-key").val();
  let endpoint = $("#ai-endpoint").val();

  if (!modelName && provider.defaultModel) {
    modelName = provider.defaultModel;
    $("#ai-model").val(modelName);
    storeModelName(providerId, modelName);
  }

  if ((providerId === "custom" && !endpoint) || (providerId !== "seqera" && !modelName) || !aiApiKey) {
    // Open the AI toolbox section, focus on missing fields
    const endpointGroup = $("#ai_endpoint_group");
    const endpointInput = endpointGroup.find("input");
    const modelGroup = $("#ai_model_group");
    const modelInput = modelGroup.find("input");
    const apiKeyGroup = $("#ai_api_key_group");
    const apiKeyInput = apiKeyGroup.find("input");

    if (providerId === "custom" && !endpoint) {
      const endpointGroup = $("#ai_endpoint_group");
      const endpointInput = endpointGroup.find("input");
      const label = endpointGroup.find("label");
      const originalLabelColor = label.css("color");
      label.css("color", "#a94442");
      endpointGroup.addClass("has-error");
      // Remove has-error class when user makes a selection
      endpointInput.one("change", function () {
        endpointGroup.removeClass("has-error");
        label.css("color", originalLabelColor);
      });
    }

    if (!modelName) {
      const label = modelGroup.find("label");
      modelGroup.addClass("has-error");
      const originalLabelColor = label.css("color");
      label.css("color", "#a94442");
      modelInput.focus();
      // Remove has-error class when user makes a selection
      modelInput.one("change", function () {
        modelGroup.removeClass("has-error");
        label.css("color", originalLabelColor);
      });
    }

    if (!aiApiKey) {
      const label = apiKeyGroup.find("label");
      apiKeyGroup.addClass("has-error");
      const required = apiKeyGroup.find("#ai_api_key_info_required");
      const originalRequiredColor = required.css("color");
      const originalLabelColor = label.css("color");
      required.css("color", "#a94442");
      label.css("color", "#a94442");
      apiKeyInput.focus();
      // Remove has-error class when user starts typing
      apiKeyInput.one("input", function () {
        apiKeyGroup.removeClass("has-error");
        required.css("color", originalRequiredColor);
        label.css("color", originalLabelColor);
      });
    }

    mqc_toolbox_openclose("#mqc_ai", true);
    if (providerId === "custom" && !endpoint) {
      endpointInput.focus();
    } else if (!modelName) {
      modelInput.focus();
    } else if (!aiApiKey) {
      apiKeyInput.focus();
    }
    return;
  }

  let budget = getPromptBudget(modelName);
  let content = buildContent(budget);
  const totalTokens = countTokens(systemPrompt + content);
  if (totalTokens > budget) {
    errorDiv
      .text(
        `The report data is too large for the context window of ${modelName} (about ${totalTokens} > ${budget} tokens)`,
      )
      .show();
    if (wrapperDiv) wrapperDiv.show();
    return;
  }

  // Disable button and show loading state
  const buttonTitle = providerId === "custom" ? "custom endpoint" : provider.name;
  button.prop("disabled", true).html(`Requesting ${buttonTitle}...`);

  function wrapUpResponse() {
    const endpoint = $("#ai-endpoint").val();
    const disclaimerProvider = providerId === "custom" ? endpoint : provider.name;
    disclaimerDiv.find(".ai-summary-disclaimer-provider").text(disclaimerProvider);
    disclaimerDiv.find(".ai-summary-disclaimer-model").text(modelName);
    disclaimerDiv.show();
    button.data("action", "clear").prop("disabled", false).html(clearText).addClass("ai-local-content");
  }

  const startTime = performance.now();
  // If the provider says the prompt is too long, rebuild it from less of the report and try again
  let overflowRetries = 3;
  const generate = () => {
    let receievedMarkdown = "";
    window.runStreamGeneration({
      title: title + `, created on ${configCreationDate}`,
      systemPrompt: systemPrompt,
      userMessage: content,
      tags: ["multiqc"],
      onStreamStart: (resolvedModelName) => {
        modelName = resolvedModelName;
        button.html(`Starting generation...`);
      },
      onStreamNewToken: (token) => {
        responseDiv.show();
        if (wrapperDiv) wrapperDiv.show();
        receievedMarkdown += token;
        responseDiv.html(window.markdownToHtml(receievedMarkdown));
        button.html(`Generating...`);
      },
      onStreamRestart: () => {
        // The model declined and another model answers from the start
        receievedMarkdown = "";
        responseDiv.html("");
      },
      onStreamError: (error) => {
        const overflow = receievedMarkdown ? null : parseContextOverflow(String(error));
        if (overflow && (isGlobal || sectionAnchor === "general_stats_table")) {
          if (overflowRetries-- > 0) {
            const estimatedTokens = countTokens(systemPrompt + content);
            budget = shrinkBudget(budget, estimatedTokens, overflow);
            content = buildContent(budget);
            button.html(`Retrying with less of the report...`);
            generate();
            return;
          }
        }
        showToast("Error generating summary", error, "error");
        wrapUpResponse(disclaimerDiv, provider.name, modelName);
        if (!isMore && isGlobal) $("#global_ai_summary_more_button_and_disclaimer").hide();
        disclaimerDiv.hide();
      },
      onStreamComplete: (threadId) => {
        wrapUpResponse(disclaimerDiv, provider.name, modelName);
        // Update the "Chat with Seqera AI" button to point to new thread
        if (threadId) {
          continueInChatButton.attr("href", `${seqeraWebsite}/chat/${threadId}`).show();
        }
        // Enable tooltips in the new content
        const tooltipTriggerList = responseDiv.find('[data-bs-toggle="tooltip"]');
        tooltipTriggerList.each(function () {
          new bootstrap.Tooltip(this);
        });

        // Save response to localStorage
        const elementId = button.data("plot-anchor") || "global";
        localStorage.setItem(
          `ai_response_${reportUuid}_${elementId}${isMore ? "_more" : ""}`,
          JSON.stringify({
            text: receievedMarkdown,
            provider: providerId,
            model: modelName,
            timestamp: Date.now(),
            threadId: threadId,
            endpoint: endpoint,
          }),
        );
        const endTime = performance.now();
        console.log(`Time to generate more: ${endTime - startTime}ms`);
        if (!isMore && isGlobal) $("#global_ai_summary_more_button_and_disclaimer").show();
        disclaimerDiv.show();
      },
    });
  };
  generate();
}

// Global (report-level) summary generation
async function generateCallback(e) {
  e.preventDefault();
  const button = $(e.currentTarget);
  const isMore = button.hasClass("ai-generate-button-more");
  const action = button.data("action");
  const responseDiv = $("#" + button.data("response-div"));
  const detailedAnalysisDiv = $("#" + button.data("detailed-analysis-div"));
  const errorDiv = $("#" + button.data("error-div"));
  const wrapperDiv = $("#" + button.data("wrapper-div"));
  const originalButtonHtml = button.data("original-html");
  const elementId = button.data("plot-anchor") || "global";
  const continueInChatButton = $("#" + button.data("continue-in-chat-button"));

  if (action === "clear") {
    e.preventDefault();
    localStorage.removeItem(`ai_response_${reportUuid}_${elementId}${isMore ? "_more" : ""}`);
    responseDiv.html("").hide();
    detailedAnalysisDiv.html("").hide();
    errorDiv.html("").hide();
    if (wrapperDiv) wrapperDiv.hide();
    continueInChatButton.hide();
    button.html(originalButtonHtml).data("action", "generate").removeClass("ai-local-content");
  } else {
    summarizeWithAi(button);
  }
}

$(function () {
  // Summaries generated at report build time are stored as markdown and rendered here,
  // so they go through exactly the same conversion and sanitising as ones generated in
  // the browser. Escaped in the HTML at rest, so textContent gives back the markdown.
  $(".ai-summary-markdown").each(function () {
    const markdown = this.textContent.trim();
    if (markdown) this.innerHTML = window.markdownToHtml(markdown);
  });

  $("#global_ai_summary_expand").each(function () {
    const responseDiv = $("#global_ai_summary_detailed_analysis_response");
    const isLocalContent = responseDiv.hasClass("ai-local-content");

    const expandBtn = $("#global_ai_summary_expand");
    const expandBtnSvg = expandBtn.find("svg");
    const expandBtnText = expandBtn.find("span");

    let isExpanded = responseDiv.is(":visible");
    const storedState = localStorage.getItem("mqc_ai_global_summary_expanded");
    if (storedState === "expanded") isExpanded = true;
    if (storedState === "collapsed") isExpanded = false;

    if (isExpanded && !isLocalContent) {
      responseDiv.show();
      expandBtnSvg.css("transform", "rotate(180deg)");
      expandBtnText.text("Hide full summary");
    } else {
      responseDiv.hide();
      expandBtnSvg.css("transform", "rotate(0deg)");
      expandBtnText.text("View full summary");
    }

    expandBtn.on("click", (e) => {
      e.preventDefault();
      isExpanded = !isExpanded;
      if (isExpanded) {
        responseDiv.show();
        expandBtnSvg.css("transform", "rotate(180deg)");
        expandBtnText.text("Hide full summary");
      } else {
        responseDiv.hide();
        expandBtnSvg.css("transform", "rotate(0deg)");
        expandBtnText.text("View full summary");
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
    const continueInChatButton = $("#" + button.data("continue-in-chat-button"));
    if (wrapperDiv) wrapperDiv.addClass("ai-local-content");

    if (action === "clear") {
      button.html(clearText).addClass("ai-local-content");
    } else {
      // Load cached AI responses on page load
      const plotAnchor = button.data("plot-anchor") || "global";
      const cachedSummaryDump = localStorage.getItem(`ai_response_${reportUuid}_${plotAnchor}${isMore ? "_more" : ""}`);
      if (cachedSummaryDump) {
        const cachedSummary = JSON.parse(cachedSummaryDump);
        responseDiv.show().html(window.markdownToHtml(cachedSummary.text));
        if (wrapperDiv) wrapperDiv.show();
        const provider = window.AI_PROVIDERS[cachedSummary.provider];
        disclaimerDiv
          .find(".ai-summary-disclaimer-provider")
          .text(provider.name == "Custom" ? cachedSummary.endpoint : provider.name);
        disclaimerDiv.find(".ai-summary-disclaimer-model").text(cachedSummary.model);
        disclaimerDiv.show();
        button.html(clearText).data("action", "clear").prop("disabled", false).addClass("ai-local-content");

        const threadId = cachedSummary.threadId;
        if (threadId) {
          continueInChatButton.attr("href", `${seqeraWebsite}/chat/${threadId}`);
          continueInChatButton.show();
        }
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
      $(".mqc_regex_mode input").prop("checked", true);
    }

    let highlightedSamples = window.mqc_highlight_f_texts;
    if (!highlightedSamples.includes(sampleName)) {
      $("#mqc_colour_filter").val(sampleName);
      // Use the current color from the palette (set by highlights.js form handler)
      $(this).css("font-weight", "bold");
      // also highlight all <sample> elements in text that match the sample name
      $("sample").each(function () {
        if ($(this).text().indexOf(sampleName) > -1) $(this).css("font-weight", "bold");
      });
      // Only submit form when adding a new highlight
      $("#mqc_color_form").trigger("submit");
    } else {
      $("#mqc_col_filters li").each(function () {
        if ($(this).find(".f_text").val() === sampleName) {
          $(this).find(".btn-close").click();
        }
      });
      $(this).css("font-weight", "normal");
      // also remove the bold from all <sample> elements in text that match the sample name
      $("sample").each(function () {
        if ($(this).text().indexOf(sampleName) > -1) $(this).css("font-weight", "normal");
      });
    }

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
    const text = window.multiqcDescription + "\n" + systemPrompt + "\n\n" + content;

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
  // Character-based estimate. Report data is mostly numbers and separators, which tokenize at
  // about 2 characters per token, rather than the ~4 of prose. If the provider says the prompt
  // is too long, the prompt is rebuilt using the provider's own count.
  return Math.ceil(text.length / 2);
}
