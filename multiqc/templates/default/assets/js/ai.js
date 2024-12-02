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

async function formatPlotForAi(button) {
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

  return reportData;
}

// Global (report-level) summary generation
async function generateCallback(e) {
  e.preventDefault();

  const button = $(e.currentTarget);
  const sectionAnchor = button.data("section-anchor") || "global";

  let content;
  let systemPrompt;
  if (sectionAnchor === "global") {
    content = button.data("content-base64");
    content = atob(content);
    systemPrompt = systemPromptReport;
  } else {
    content = await formatPlotForAi(button);
    systemPrompt = systemPromptPlot;
  }

  responseDiv = $("#" + sectionAnchor + "_ai_detailed_summary");
  errorDiv = $("#" + sectionAnchor + "_ai_summary_error");
  disclaimerDiv = $("#" + sectionAnchor + "_ai_summary_disclaimer");
  wrapperDiv = $("#" + sectionAnchor + "_ai_summary");

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
