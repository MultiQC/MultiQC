////////////////////////////////////////////////
// MultiQC Report Toolbox - AI Functionality
////////////////////////////////////////////////

function updatePanel(providerId) {
  const provider = AI_PROVIDERS[providerId];

  // Add label dynamically
  let aiModelInfo = "";
  let aiApiKeyInfo = "";

  // Update model field with stored or default value for new provider
  if (providerId === "seqera") {
    $(".ai-generate-button-wrapper").show();
    $(".ai-copy-button-wrapper").hide();
    $("#ai_provider_info").html("");
    $("#ai_model_group").hide();
    $("#ai_endpoint_group").hide();
    $("#ai_query_options_group").hide();
    $("#ai_api_key_group").show();
    $("#ai_api_key_info").html(
      `<span id="ai_api_key_info_required">Field is required.</span> You can find your API key in the <a style='text-decoration: underline;' href='${provider.apiKeysUrl}' target='_blank'>${provider.name} dashboard</a>. `,
    );
    $(".ai_provider_logos").hide();
    $("#ai_provider_logo_seqera").show();
    $("#mqc_anonymize_samples_switch").show();
    $("#ai_context_window_group").hide();
  } else if (providerId === "clipboard") {
    $(".ai-generate-button-wrapper").hide();
    $(".ai-copy-button-wrapper").show();
    $("#ai_model_group").hide();
    $("#ai_api_key_group").hide();
    $("#ai_provider_info").html(
      `Copy AI LLM prompt, including report data, to the clipboard. For use with any 3rd party AI tools.`,
    );
    $("#ai_endpoint_group").hide();
    $("#ai_query_options_group").hide();
    $(".ai_provider_logos").hide();
    $("#mqc_anonymize_samples_switch").show();
    $("#ai_context_window_group").hide();
  } else if (providerId === "none") {
    $(".ai-generate-button-wrapper").hide();
    $(".ai-copy-button-wrapper").hide();
    $("#ai_provider_info").html("Remove 'Summarize' buttons from report if you're not going to use AI summaries.");
    $(".ai_provider_logos").hide();
    $("#ai_model_group").hide();
    $("#ai_api_key_group").hide();
    $("#ai_endpoint_group").hide();
    $("#ai_query_options_group").hide();
    $("#mqc_anonymize_samples_switch").hide();
    $("#ai_context_window_group").hide();
  } else {
    // OpenAI, Anthropic, custom - showing logos and more helpful guides
    $(".ai-generate-button-wrapper").show();
    $(".ai-copy-button-wrapper").hide();
    $("#ai_provider_info").html("");
    $(".ai_provider_logos").hide();
    $("#ai_model_group").show();
    $("#ai_api_key_group").show();
    $("#mqc_anonymize_samples_switch").show();
    if (providerId === "custom") {
      $("#ai_endpoint_group").show();
      $("#ai_query_options_group").show();
      $("#ai_context_window_group").show();
      aiApiKeyInfo = `<span id="ai_api_key_info_required">Field is required.</span></a>`;
      $("#ai_model_info").hide();
      $("#ai_api_key_info").hide();
    } else {
      $("#ai_endpoint_group").hide();
      $("#ai_query_options_group").hide();
      $("#ai_context_window_group").hide();
      aiApiKeyInfo = `<span id="ai_api_key_info_required">Field is required.</span> You can find your API key in the <a style='text-decoration: underline;' href='${provider.apiKeysUrl}' target='_blank'>${provider.name} console</a>. `;
      // Add clickable model suggestions if available
      let suggestedModels = provider.suggestedModels || [];
      aiModelInfo = `You can find the available models for ${provider.name} in the <a style='text-decoration: underline;' href='${provider.modelsUrl}' target='_blank'>${provider.name} documentation</a>.`;
      if (suggestedModels.length > 0) {
        aiModelInfo += " For example: ";
        aiModelInfo += suggestedModels
          .map((model) => `<a href="#" class="ai-model-suggestion" data-model="${model}"><code>${model}</code></a>`)
          .join(", ");
      }
      $("#ai_model_info").html(aiModelInfo).show();
      $("#ai_api_key_info").html(aiApiKeyInfo).show();
    }
    // Doing it here again because model depends on provider
    const storedModel = getStoredModelName(providerId);
    const defaultModel = provider && provider.defaultModel ? provider.defaultModel : null;
    $("#ai-model").val(storedModel || defaultModel);

    $(".ai_provider_logos").hide();
    if (providerId === "openai") {
      $("#ai_provider_logo_openai").show();
    } else if (providerId === "anthropic") {
      $("#ai_provider_logo_anthropic").show();
    }
  }
}

// Storage functions
function getStoredProvider() {
  let storedProviderId = getFromLocalStorage("mqc_ai_provider");
  if (storedProviderId && window.AI_PROVIDERS[storedProviderId]) return storedProviderId;
  return null;
}
function storeProvider(providerId) {
  saveToLocalStorage("mqc_ai_provider", providerId);
}
function getStoredApiKey(providerId) {
  return getFromLocalStorage(`mqc_ai_key_${providerId}`);
}
function storeApiKey(providerId, apiKey) {
  saveToLocalStorage(`mqc_ai_key_${providerId}`, apiKey);
}
function getStoredModelName(providerId) {
  return getFromLocalStorage(`mqc_ai_model_${providerId}`);
}
function storeModelName(providerId, modelName) {
  saveToLocalStorage(`mqc_ai_model_${providerId}`, modelName);
}
function getStoredEndpoint() {
  return getFromLocalStorage(`mqc_ai_endpoint`);
}
function storeEndpoint(endpoint) {
  saveToLocalStorage(`mqc_ai_endpoint`, endpoint);
}
function getStoredQueryOptions() {
  let data = getFromLocalStorage(`mqc_ai_query_options`);
  if (data) {
    try {
      return JSON.parse(data);
    } catch (e) {
      console.error("Error parsing extra query options", e);
      return null;
    }
  }
  return null;
}
function storeQueryOptions(options) {
  saveToLocalStorage(`mqc_ai_query_options`, JSON.stringify(options));
}
function getStoredSampleAnonymizationEnabled() {
  const stored = getFromLocalStorage("mqc_ai_anonymize_samples");
  return stored === null ? null : stored === "true";
}
function storeSampleAnonymizationEnabled(value) {
  saveToLocalStorage("mqc_ai_anonymize_samples", value.toString());
}
function getStoredContextWindow() {
  return getFromLocalStorage(`mqc_ai_context_window`);
}

// Make function available globally
window.getStoredSampleAnonymizationEnabled = getStoredSampleAnonymizationEnabled;

// Make functions available globally
window.initAI = function () {
  // Populate provider dropdown dynamically
  const aiProviderSelect = $("#ai-provider");
  aiProviderSelect.empty();
  Object.entries(AI_PROVIDER_GROUPS).forEach(([groupName, groupProviders]) => {
    let optgroup = $("<optgroup>", { label: groupName });
    groupProviders.forEach((groupProviderID) => {
      optgroup.append(
        $("<option>", {
          value: groupProviderID,
          text: window.AI_PROVIDERS[groupProviderID].name,
        }),
      );
    });
    aiProviderSelect.append(optgroup);
  });

  // Set initial values from storage or values from Python
  const providerId = getStoredProvider() || aiConfigProviderId || "seqera";
  aiProviderSelect.val(providerId);
  const provider = window.AI_PROVIDERS[providerId];
  $("#ai-api-key").val(getStoredApiKey(providerId) || "");

  let model = getStoredModelName(providerId);
  if (model === null && aiConfigModel !== "None") model = aiConfigModel;
  if (model === null && provider.defaultModel) model = provider.defaultModel;
  $("#ai-model").val(model);

  let endpoint = getStoredEndpoint();
  if (endpoint === null) endpoint = aiConfigCustomEndpoint;
  $("#ai-endpoint").val(endpoint);

  let storedOptions = getStoredQueryOptions();
  if (storedOptions === null) storedOptions = aiConfigExtraQueryOptions;
  $("#ai-query-options").val(JSON.stringify(storedOptions));

  let storedContextWindow = getStoredContextWindow();
  if (storedContextWindow === null && aiConfigCustomContextWindow !== "None") {
    try {
      storedContextWindow = parseInt(aiConfigCustomContextWindow);
    } catch {}
  }
  storedContextWindow = storedContextWindow ?? 128000;
  $("#ai-context-window").val(storedContextWindow);

  updatePanel(providerId);

  // Save provider changes
  aiProviderSelect.change(function () {
    const selectedProviderId = $(this).val();
    storeProvider(selectedProviderId);

    // Update API key field
    const storedKey = getStoredApiKey(selectedProviderId);
    $("#ai-api-key").val(storedKey || "");

    updatePanel(selectedProviderId);
  });

  // Save model changes
  $("#ai-model").change(function () {
    const providerId = $("#ai-provider").val();
    const model = $(this).val();
    storeModelName(providerId, model);
  });

  $("#ai-endpoint").change(function () {
    const endpoint = $(this).val();
    storeEndpoint(endpoint);
  });

  // Save API key changes
  $("#ai-api-key").change(function () {
    const providerId = $("#ai-provider").val();
    const apiKey = $(this).val();
    storeApiKey(providerId, apiKey);
  });

  // Initialize anonymize samples switch
  const anonymizeSamplesEnabled = getStoredSampleAnonymizationEnabled() ?? aiAnonymizeSamples === "true";
  if (anonymizeSamplesEnabled) {
    $(".mqc_switch.anonymize_samples").removeClass("off").addClass("on").text("on");
  }

  // Handle anonymize samples switch clicks
  $("#mqc_anonymize_samples").on("change", function () {
    storeSampleAnonymizationEnabled($(this).prop("checked"));
  });
};

// Initialize AI cookie reading
window.initAICookies = function () {
  // Read the JWT local cookie in in the seqera.io domain - that's our API key:
  const jwt = document.cookie.split("; ").find((row) => row.startsWith("jwt="));
  if (jwt) {
    const jwtValue = jwt.split("=")[1];
    saveToLocalStorage(`mqc_ai_key_seqera`, jwtValue);
  }
};
