////////////////////////////////////////////////
// AI stuff
////////////////////////////////////////////////

window.continueInChatHandler = function (event) {
  let el = $(event.currentTarget);
  let seqeraWebsite = el.data("seqera-website");

  // Either report uuid, or encoded system and chat messages
  let generationId = el.data("generation-id");
  let encodedSystemMessage = el.data("encoded-system-message");
  let encodedChatMessages = el.data("encoded-chat-messages");

  let url = seqeraWebsite + "/ask-ai/";
  if (generationId) {
    url += "?generation-id=" + generationId;
  }

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

async function globalGenerateButtonCallback(e) {
  e.preventDefault();

  const button = $(e.currentTarget);
  const contentBase64 = button.data("content-base64");
  const content = atob(contentBase64);

  const wrapperDiv = $("#global_ai_summary");
  const responseDiv = $("#global_ai_detailed_summary");
  const errorDiv = $("#global_ai_summary_error");
  const disclaimerDiv = $("#global_ai_summary_disclaimer");

  generateWithLLM(
    button,
    responseDiv,
    errorDiv,
    disclaimerDiv,
    wrapperDiv,
    content,
    systemPromptReport,
    globalGenerateButtonCallback,
    null,
    "Generate more details...",
  );
}

async function plotGenerateSummaryReset(button) {
  let sectionAnchor = button.data("section-anchor");
  const wrapperDiv = $("#" + sectionAnchor + "_ai_summary");
  wrapperDiv.hide();
}

async function plotGenerateSummaryButtonCallback(e) {
  e.preventDefault();
  const button = $(e.currentTarget);

  let moduleAnchor = button.data("module-anchor");
  let sectionAnchor = button.data("section-anchor");
  let plotAnchor = button.data("plot-anchor");

  let plot = mqc_plots[plotAnchor];
  let formattedData = plot.prepDataForLlm();

  let moduleWrapper = $("#mqc-module-section-" + moduleAnchor);
  let moduleName = moduleWrapper.find(".mqc-module-title").text();
  let moduleSummary = moduleWrapper.find(".mqc-module-section-first > p").text();
  let moduleComment = moduleWrapper.find(".mqc-section-comment").text();

  let sectionName = button.data("section-name");
  let sectionWrapper = $("#mqc-section-wrapper-" + sectionAnchor);
  let description = sectionWrapper.find(".mqc-section-description").text();
  let comment = sectionWrapper.find(".mqc-section-comment").text();
  let helptext = sectionWrapper.find(".mqc-section-helptext").text();

  // Elements to populate
  const responseDiv = $("#" + sectionAnchor + "_ai_detailed_summary");
  const disclaimerDiv = $("#" + sectionAnchor + "_ai_summary_disclaimer");
  const errorDiv = $("#" + sectionAnchor + "_ai_summary_error");
  const wrapperDiv = $("#" + sectionAnchor + "_ai_summary");

  let reportData = `
QC tool that generated data for the plot: ${moduleName} (${moduleSummary}). ${
    moduleComment ? `\nComment: ${moduleComment}` : ""
  }

Section: ${sectionName} ${description ? `(${description})` : ""}. ${comment ? `\nComment: ${comment}` : ""} ${
    helptext ? `\nHelptext: ${helptext}` : ""
  }

Plot type: ${plot.plotType}. ${plot.plotDescription ? `\nPlot description: ${plot.plotDescription}` : ""} ${
    plot.plotHelptext ? `\nPlot helptext: ${plot.plotHelptext}` : ""
  }

Plot data:
${formattedData}
`;

  generateWithLLM(
    button,
    responseDiv,
    errorDiv,
    disclaimerDiv,
    wrapperDiv,
    reportData,
    systemPromptPlot,
    plotGenerateSummaryButtonCallback,
    plotGenerateSummaryReset,
    "AI summary",
  );
}

$(function () {
  // "Show More" button to expand pre-generated full AI summary
  $(".ai-summary").each(function () {
    const $details = $(this).find("details");
    const $showMoreBtn = $(this).find(".ai-summary-expand");
    $showMoreBtn.on("click", function (e) {
      if ($details.prop("open")) {
        $details.prop("open", false);
      } else {
        $details.prop("open", true);
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
    let elementId = null;
    let generateCallback = null;
    let resetCallback = null;
    let responseDiv = null;
    let disclaimerDiv = null;
    let errorDiv = null;
    let wrapperDiv = null;
    if (button.data("plot-anchor")) {
      sectionAnchor = button.data("section-anchor");
      elementId = button.data("plot-anchor");
      generateCallback = plotGenerateSummaryButtonCallback;
      resetCallback = plotGenerateSummaryReset;
      responseDiv = $("#" + sectionAnchor + "_ai_detailed_summary").show();
      disclaimerDiv = $("#" + sectionAnchor + "_ai_summary_disclaimer").show();
      errorDiv = $("#" + sectionAnchor + "_ai_summary_error");
      wrapperDiv = $("#" + sectionAnchor + "_ai_summary");
    } else {
      elementId = "global";
      generateCallback = globalGenerateButtonCallback;
      resetCallback = null;
      responseDiv = $("#global_ai_detailed_summary");
      disclaimerDiv = $("#global_ai_summary_disclaimer");
      errorDiv = $("#global_ai_summary_error");
      wrapperDiv = $("#global_ai_summary");
    }
    const originalButtonText = button.text();
    const cachedResponse = localStorage.getItem(`ai_response_${reportUuid}_${elementId}`);
    if (cachedResponse) {
      // Load cached AI responses on page load
      const responseData = JSON.parse(cachedResponse);
      responseDiv.html(markdownToHtml(responseData.text));
      if (wrapperDiv) wrapperDiv.show();
      wrapUpResponse(
        responseDiv,
        disclaimerDiv,
        errorDiv,
        wrapperDiv,
        button,
        responseData.provider,
        responseData.model,
        generateCallback,
        resetCallback,
        originalButtonText,
      );
    } else {
      button.click(generateCallback);
    }
  });

  // Click handler to highlight samples
  $(document).on("click", "sample", function (e) {
    e.preventDefault();
    let sampleName = $(this).text();
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
});

async function wrapUpResponse(
  responseDiv,
  disclaimerDiv,
  errorDiv,
  wrapperDiv,
  button,
  provider,
  model,
  generateCallback,
  resetCallback,
  originalButtonText,
) {
  disclaimerDiv.html(`This summary is AI-generated. Provider: ${provider}, model: ${model}`).show();
  const elementId = button.data("plot-anchor") || "global";
  button
    .text("Reset local content")
    .prop("sa", "background-color: #f2f2f2;")
    .prop("disabled", false)
    .off("click")
    .on("click", function (e) {
      e.preventDefault();
      if (resetCallback) resetCallback(button);
      localStorage.removeItem(`ai_response_${reportUuid}_${elementId}`);
      responseDiv.html("");
      errorDiv.html("");
      disclaimerDiv.html("");
      wrapperDiv.hide();
      button.text(originalButtonText).prop("style", "background-color: white;").off("click").click(generateCallback);
    });
}

async function generateWithLLM(
  button,
  responseDiv,
  errorDiv,
  disclaimerDiv,
  wrapperDiv,
  content,
  systemPrompt,
  generateCallback,
  resetCallback,
  originalButtonText,
) {
  // Check for stored API key
  let provider = getStoredProvider();
  let aiApiKey = getStoredApiKey(provider);
  if (!aiApiKey || aiApiKey === undefined) {
    // Open the AI toolbox section
    mqc_toolbox_openclose("#mqc_ai", true);
    return;
  }

  // Disable button and show loading state
  button.prop("disabled", true);
  button.text("Generating...");

  const startTime = performance.now();
  let fullModelName = null;
  await generateDetailedSummary();
  async function generateDetailedSummary() {
    let receievedMarkdown = "";
    streamGeneration(
      function onStreamStart(model) {
        fullModelName = model;
        wrapperDiv.show();
      },
      function onStreamNewToken(token) {
        receievedMarkdown += token;
        responseDiv.html(markdownToHtml(receievedMarkdown));
      },
      function onStreamError(error) {
        errorDiv.html(error).show();
        wrapperDiv.show();
      },
      function onStreamComplete() {
        const provider = getStoredProvider();
        wrapUpResponse(
          responseDiv,
          disclaimerDiv,
          errorDiv,
          wrapperDiv,
          button,
          provider,
          fullModelName,
          generateCallback,
          resetCallback,
          originalButtonText,
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
      systemPrompt,
      content,
      ["multiqc"],
    );
  }
}
