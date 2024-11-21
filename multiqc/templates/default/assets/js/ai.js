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

$(function () {
  // Add settings button and modal
  $("body").append(getSettingsModalHtml());

  // Set initial values when opening modal
  const provider = getStoredProvider() || "Anthropic";
  $("#ai-provider").val(provider);
  const model = getStoredModelName(provider) || AI_PROVIDERS[provider].defaultModel;
  $("#ai-model").val(model);
  const apiKey = getStoredApiKey(provider);
  $("#ai-api-key").val(apiKey || "");

  // Handle settings button click
  $(".ai-settings-button").click(function (e) {
    e.preventDefault();
    $("#aiSettingsModal").modal("show");
  });

  // Add handler to update model when provider changes
  $("#ai-provider").change(function () {
    const provider = $(this).val();
    const defaultModel = AI_PROVIDERS[provider].defaultModel;
    $("#ai-model").val(defaultModel);
  });

  // Handle settings save
  $("#saveAiSettings").click(function () {
    // Save API key for selected provider
    const provider = $("#ai-provider").val();
    storeProvider(provider);
    const model = $("#ai-model").val();
    storeModelName(provider, model);
    const apiKey = $("#ai-api-key").val();
    storeApiKey(provider, apiKey);
    $("#aiSettingsModal").modal("hide");
  });

  // Add handler for Enter key in the modal
  $("#aiSettingsModal").on("keypress", function (e) {
    if (e.which === 13) {
      $("#saveAiSettings").click();
    }
  });

  // // Add this to update the API key field when provider changes
  // $("#defaultProvider").change(function () {
  //   const provider = $(this).val();
  //   const storedKey = getStoredApiKey(provider);
  //   $("#ai-token").val(storedKey || "");
  // });

  // Add "Show More" button to AI summary
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

  $(".ai-summary sample").hover(function () {
    $(this).css("opacity", 0.9);
  });

  $(".ai-generate-more-plot").click(async function (e) {
    e.preventDefault();
    let button = $(this);
    let buttonContainer = button.parent();
    // Disable button and show spinner
    button.prop("disabled", true);
    button.html("Generating...");

    let moduleAnchor = button.data("module-anchor");
    let sectionAnchor = button.data("section-anchor");
    let plotAnchor = button.data("plot-anchor");
    let multiqcVersion = button.data("multiqc-version");
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

    const responseDiv = $("#" + sectionAnchor + "_ai_summary .ai-summary-main-text");

    const errorDiv = $("#" + sectionAnchor + "_ai_summary .ai-summary-error");

    let fullModelName = null;

    const startTime = performance.now();
    await generateDetailedSummary();
    async function generateDetailedSummary() {
      let receievedText = "";
      streamGeneration(
        function onStreamStart(model) {
          $("#" + sectionAnchor + "_ai_summary").show();
          fullModelName = model;
        },
        function onStreamNewToken(token) {
          receievedText += token;
          responseDiv.html(markdownToHtml(receievedText));
        },
        function onStreamError(error) {
          errorDiv.html(error);
        },
        function onStreamComplete() {
          const provider = getStoredProvider();
          responseDiv.append(
            `<div class="ai-summary-disclaimer">This summary is AI-generated. 
            Provider: ${provider}, model: ${fullModelName}</div>`,
          );
          const endTime = performance.now();
          console.log(`Time to generate more: ${endTime - startTime}ms`);
          buttonContainer.hide();
        },
        systemPromptPlot,
        reportData,
        ["multiqc", "multiqc-plot-summary", `multiqc-${multiqcVersion}`],
      );
    }
  });

  // Add click event listeners to clickable samples
  $(document).on("click", "sample", function (e) {
    e.preventDefault();
    var sampleName = $(this).text();
    var color = $(this).css("color");
    var highlightedSamples = window.mqc_highlight_f_texts;
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

  // Add this to the $(function() {...}) block to set the initial selected provider
  $("#defaultProvider").val(sessionStorage.getItem("multiqc_default_ai_provider") || "Anthropic");

  // Handle modal shown event
  $("#aiSettingsModal").on("shown.bs.modal", function () {
    // Focus first empty input field
    const inputs = $("#aiSettingsModal").find("input");
    const firstEmptyInput = inputs
      .filter(function () {
        return !this.value;
      })
      .first();

    if (firstEmptyInput.length) {
      firstEmptyInput.focus();
    } else {
      inputs.first().focus();
    }
  });
});

async function generateMoreHandler(event) {
  event.preventDefault();
  let button = $(event.currentTarget);

  // Check for stored API key
  let provider = getStoredProvider();
  let aiApiKey = getStoredApiKey(provider);
  if (!aiApiKey || aiApiKey === undefined) {
    // Add one-time handler for the save button
    const saveHandler = function () {
      const provider = $("#ai-provider").val();
      storeProvider(provider);
      const apiKey = $("#ai-api-key").val();
      if (apiKey) {
        storeApiKey(provider, apiKey);
        $("#aiSettingsModal").modal("hide");
        // Remove this one-time handler
        $("#saveAiSettings").off("click", saveHandler);

        // Continue with generation
        generateMore();
      }
    };
    $("#saveAiSettings").on("click", saveHandler);
    $("#aiSettingsModal").modal("show");
    return;
  } else {
    generateMore();
  }

  async function generateMore() {
    // Disable button and show loading state
    button.prop("disabled", true);
    button.text("Generating...");
    let buttonContainer = button.parent();

    const contentBase64 = button.data("content-base64");
    const content = atob(contentBase64);

    const summaryDiv = buttonContainer.prev(".ai-short-summary");

    let fullModelName = null;

    const startTime = performance.now();
    await generateDetailedSummary();
    async function generateDetailedSummary() {
      const responseDiv = $('<div class="ai-detailed-summary"></div>');
      summaryDiv.after(responseDiv);

      const errorDiv = $('<div class="ai-summary-error"></div>');
      responseDiv.after(errorDiv);

      let receievedText = "";
      streamGeneration(
        function onStreamStart(model) {
          buttonContainer.hide();
          fullModelName = model;
        },
        function onStreamNewToken(token) {
          receievedText += token;
          responseDiv.html(markdownToHtml(receievedText));
        },
        function onStreamError(error) {
          errorDiv.html(error);
        },
        function onStreamComplete() {
          const provider = getStoredProvider();
          responseDiv.append(
            `<div class="ai-summary-disclaimer">This summary is AI-generated. 
            Provider: ${provider}, model: ${fullModelName}</div>`,
          );
          const endTime = performance.now();
          console.log(`Time to generate more: ${endTime - startTime}ms`);
        },
        systemPromptReport,
        content,
        ["multiqc", "multiqc-generate-more"],
      );
    }
  }
}
