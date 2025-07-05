////////////////////////////////////////////////
// MultiQC Report Toolbox - AI Functionality
////////////////////////////////////////////////

function updatePanel(providerId) {
  const provider = window.AI_PROVIDERS[providerId];

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
      `<span id="ai_api_key_info_required">Field is required.</span> You can find your API key in the <a style='text-decoration: underline;' href='${provider.apiKeysUrl}' target='_blank'>${provider.name} dashboard</a>`,
    );
    $("#ai_provider_logo").html(
      '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 386.21 80" height="22" width="112" class="_logo_1pric_80"><defs><clipPath id="Logo_svg__a"><path d="M0 0h386.21v80H0z" style="fill:none;stroke-width:0"></path></clipPath><clipPath id="Logo_svg__b"><path d="M0 0h79.68v80H0z" style="fill:none;stroke-width:0"></path></clipPath></defs><g style="clip-path:url(#Logo_svg__a)"><g style="clip-path:url(#Logo_svg__b)"><path d="m67.65 41.42.11.04.11.04.11.05.11.06.1.06.1.06.1.07.09.07.09.08.09.08.08.09.08.09.08.1.07.1.07.1.06.1.05.11.05.11.04.11.04.11.03.12.03.11.02.12.02.11v.12l.02.12v.36l-.02.12-.02.12-.02.12-.03.12-.03.11-.04.12-.04.11-.05.11-.05.11-.06.1-.07.1-.07.1-.08.09-.08.09-.08.08-.09.09-.09.08-.09.07-.1.07-.1.07-.1.06-.11.06-.11.05-.11.05-.11.04-.11.04-.11.03-.11.03-.12.02h-.11l-.12.03h-.02.12l.22-.01h.22l.22-.03.21-.03.22-.03.21-.04.21-.04.21-.05.21-.05.21-.06.21-.07.2-.07.2-.08.2-.08.2-.09.2-.09.19-.1.19-.1.19-.11.18-.12.18-.12.18-.12.18-.13.17-.13.17-.13.17-.14.16-.15.16-.15.16-.16.15-.16.14-.16.14-.17.14-.17.13-.18.13-.18.12-.18.12-.19.11-.19.11-.19.11-.19.1-.19.09-.2.09-.2.08-.2.08-.2.07-.21.06-.21.06-.21.05-.21.05-.21.04-.22.04-.21.03-.21.02-.22.02-.22.02-.22v-.88l-.02-.22-.02-.22-.02-.21-.03-.22-.04-.22-.04-.21-.05-.21-.04-.16H48.51l19.02 5.25.13.04Z" style="fill:#160f26;stroke-width:0"></path><path d="m79.67 43.32-.02-.28-.02-.27-.03-.28-.03-.27-.04-.27-.04-.27-.05-.27-.05-.27-.06-.27-.07-.27-.07-.27-.08-.26-.08-.26-.09-.26-.1-.26-.1-.25-.11-.26-.11-.25-.11-.25-.12-.25-.13-.24-.13-.24-.14-.24-.14-.24-.14-.23-.15-.23-.15-.22-.16-.22-.05-.07h-2.38l.04.16.05.21.04.21.04.22.03.22.02.21.02.22.02.22v.88l-.02.22-.02.22-.02.22-.03.21-.04.21-.04.22-.05.21-.05.21-.06.21-.06.21-.07.21-.08.2-.08.2-.09.2-.09.2-.1.19-.11.19-.11.19-.11.19-.12.19-.12.18-.13.18-.13.18-.14.17-.14.17-.14.16-.15.16-.16.16-.16.15-.16.15-.17.14-.17.13-.17.13-.18.13-.18.12-.18.12-.18.12-.19.11-.19.1-.19.1-.2.09-.2.09-.2.08-.2.08-.2.07-.21.07-.21.06-.21.05-.21.05-.21.04-.21.04-.22.03-.21.03-.22.02h-.22l-.22.02H0v10.32h67.38l.27-.02.27-.02.27-.03.27-.03.27-.04.27-.04.27-.05.27-.06.27-.06.26-.06.26-.07.26-.08.26-.08.26-.09.26-.09.26-.1.25-.1.25-.11.25-.12.25-.12.24-.12.24-.13.24-.13.24-.14.23-.14.23-.15.23-.16.22-.16.22-.17.21-.17.21-.18.21-.18.2-.18.2-.19.19-.19.19-.2.18-.2.18-.2.18-.21.17-.21.17-.21.17-.22.16-.22.15-.22.15-.23.14-.23.14-.24.14-.24.13-.24.13-.24.12-.25.11-.25.11-.25.11-.25.1-.25.1-.26.09-.26.08-.26.08-.27.07-.26.07-.27.06-.27.05-.27.05-.27.04-.27.04-.27.03-.27.03-.27.02-.28.02-.27v-.27l.01-.28v-.28l-.01-.27Z" style="fill:#fa6863;stroke-width:0"></path><path d="m12.04 38.57-.11-.04-.11-.04-.11-.05-.11-.06-.1-.06-.1-.06-.1-.07-.09-.07-.09-.08-.09-.08-.08-.09-.08-.09-.08-.1-.07-.1-.07-.1-.06-.1-.05-.11-.05-.11-.04-.11-.04-.11-.03-.12-.03-.11-.03-.12-.02-.11v-.12l-.02-.12v-.36l.02-.12.02-.12.03-.12.03-.12.03-.11.04-.12.04-.11.05-.11.05-.1.06-.11.07-.1.07-.1.08-.09.08-.09.08-.09.09-.08.09-.08.1-.07.1-.07.1-.07.11-.06.1-.06.11-.05.11-.04.11-.04.11-.04.11-.03.12-.02.11-.02h.12l.12-.03h.04-.34l-.22.03-.22.02-.21.03-.22.03-.22.04-.21.04-.21.05-.21.05-.21.06-.21.07-.2.07-.21.07-.2.08-.2.09-.2.09-.19.1-.19.1-.18.11-.18.12-.18.12-.18.12-.18.13-.18.13-.17.14-.17.14-.16.15-.16.15-.16.15-.15.16-.14.16-.14.17-.13.17-.14.17-.12.18-.12.18-.12.18-.11.19-.11.19-.11.2-.1.19-.09.2-.09.2-.08.2-.08.2-.07.21-.06.21-.06.21-.05.21-.05.21-.04.21-.04.22-.03.21-.03.22-.02.21-.02.22v.88l.02.22.02.22.02.22.03.22.04.21.04.21.05.21.04.17h26.23l-19.02-5.25-.13-.04Z" style="fill:#160f26;stroke-width:0"></path><path d="M.01 36.67v.28l.04.27.03.27.03.27.04.27.04.27.05.27.05.27.06.27.07.27.07.27.08.26.08.26.09.26.1.26.1.25.11.26.11.25.11.25.12.24.13.24.13.24.14.24.14.23.14.23.15.23.16.23.16.22.05.07h2.38l-.04-.17-.05-.21-.04-.21-.04-.21-.03-.22-.02-.22-.02-.22-.02-.22v-.88l.02-.22.02-.21.02-.22.03-.21.04-.22.04-.21.05-.21.05-.21.06-.21.06-.21.07-.21.08-.2.08-.2.09-.2.09-.2.1-.19.11-.2.11-.19.11-.19.12-.18.12-.18.12-.18.14-.17.13-.17.14-.17.14-.16.15-.16.15-.16.15-.15.16-.14.16-.14.17-.13.17-.14.17-.12.18-.12.18-.12.18-.11.19-.11.19-.11.2-.1.19-.09.2-.09.2-.08.2-.08.2-.07.21-.06.21-.06.21-.05.21-.05.21-.04.21-.04.22-.03.21-.03.22-.02.22-.02h.22l.22-.01H79.7V23.23H12.31l-.27.02-.27.02-.27.03-.27.03-.27.04-.27.04-.27.05-.27.05-.26.06-.26.07-.26.07-.26.08-.26.08-.26.09-.26.09-.26.1-.25.1-.25.11-.25.12-.25.12-.24.12-.24.13-.24.13-.24.14-.23.15-.23.15-.22.16-.22.16-.22.16-.21.17-.21.18-.21.18-.2.18-.2.19-.19.19-.19.2-.19.2-.18.2-.18.2-.18.2-.19.2-.19.19-.19.19-.2.18-.2.18-.21.18-.21.17-.21.17-.22.16-.22.16-.22.15-.23.14-.23.14-.24.14-.24.13-.24.13-.24.12-.25.11-.25.11-.25.11-.25.1-.25.1-.26.09-.26.08-.26.08-.26.07-.26.07-.27.06-.27.05-.27.05-.27.04-.27.04-.27.03-.27.03-.27.02-.28v-.27l.03-.27v.82Z" style="fill:#f18046;stroke-width:0"></path><path d="m12.04 15.35-.11-.04-.11-.05-.11-.05-.1-.06-.11-.06-.1-.06-.1-.07-.1-.07-.09-.08-.09-.08-.08-.09-.08-.09-.08-.1-.07-.1-.06-.1-.06-.1-.05-.11-.05-.11-.04-.11-.04-.12-.04-.11-.03-.11-.03-.12-.02-.12-.02-.12v-.48l.02-.12.02-.12.03-.12.03-.12.04-.12.04-.11.04-.11.05-.11.05-.11.06-.1.06-.1.07-.1.08-.09.08-.09.08-.09.09-.08.09-.08.1-.07.1-.07.1-.07.11-.06.1-.06.11-.05.11-.04.11-.04.11-.04.11-.03.12-.02.11-.02h.12l.12-.02h.04-.36l-.22.03-.22.02-.21.03-.22.03-.22.04-.21.04-.21.05-.21.05-.21.06-.21.07-.2.07-.21.08-.2.08-.2.09-.2.09-.19.1-.19.1-.18.11-.19.12-.18.12-.18.12-.18.13-.18.13-.17.14-.17.14-.16.15-.16.15-.16.15-.15.16-.14.16-.14.17-.14.18-.13.18-.13.18-.12.18-.12.18-.12.18-.11.19-.1.19-.1.2-.09.2-.09.2-.08.21-.07.2-.07.21-.07.21-.06.21-.05.21-.05.21-.04.22-.04.22-.03.21-.03.22-.02.22-.02h67.27V0H12.31l-.27.03-.27.02-.27.02-.27.03-.27.04-.27.04-.27.05-.27.05-.26.06-.26.06-.26.07-.26.08-.26.08-.26.09-.26.09-.26.1-.25.1-.25.11-.25.12-.25.12-.24.13-.24.13-.24.13-.24.14-.23.15-.23.15-.22.16-.22.16-.22.17-.21.17-.21.18-.21.18-.2.18-.2.19-.19.19-.19.2-.19.2-.18.2-.18.21-.17.21-.17.21-.17.22-.16.22-.16.22-.15.23-.14.23-.14.23-.13.24-.13.24-.13.24-.12.25-.11.25-.11.25-.11.25-.1.25-.1.26-.09.26-.09.26-.08.26-.07.27-.07.27-.06.27-.05.27-.05.27-.04.27-.04.27-.03.27-.03.27-.02.28v.27l-.03.27v.82Z" style="fill:#160f26;stroke-width:0"></path><path d="M115.57 65.34H95.13v-8.07h20.44c3.17 0 5.6-.61 7.32-1.85 1.71-1.23 2.57-2.57 2.57-4.77s-.81-3.37-2.42-4.28c-1.62-.91-4.01-1.65-7.17-2.24l-3.3-.58c-3.23-.58-6.17-1.43-8.82-2.53s-4.75-2.63-6.3-4.57c-1.55-1.95-2.33-4.44-2.33-7.49 0-4.54 1.68-8.06 5.04-10.55 3.36-2.5 7.82-3.74 13.37-3.74h19.96v8.17h-19.96c-2.71 0-4.84.5-6.39 1.51-1.55 1-2.33 2.42-2.33 4.23 0 1.95.76 3.37 2.28 4.28q2.28 1.365 6.15 2.04l3.39.58c3.42.58 6.56 1.4 9.4 2.43 2.84 1.04 5.09 2.53 6.74 4.47 1.65 1.95 2.47 4.54 2.47 7.78 0 4.8-1.78 8.53-5.33 11.19q-5.325 3.99-14.34 3.99M238.44 79.99v-22.5h-1.55c-.71 1.3-1.74 2.55-3.1 3.74-1.36 1.2-3.09 2.19-5.18 2.97q-3.15 1.17-7.8 1.17c-4.01 0-7.69-.97-11.05-2.92s-6.04-4.75-8.04-8.42c-2-3.66-3-8.09-3-13.28v-1.46c0-5.19 1.02-9.61 3.05-13.28 2.03-3.66 4.73-6.47 8.09-8.42s7.01-2.92 10.95-2.92c4.65 0 8.22.84 10.71 2.53s4.34 3.6 5.57 5.74h1.55v-8.27h9.79v65.3h-9.98Zm-14.82-23.38c4.39 0 7.98-1.4 10.75-4.18s4.17-6.78 4.17-11.96v-.88c0-5.12-1.41-9.08-4.21-11.87-2.81-2.79-6.38-4.18-10.71-4.18s-7.8 1.4-10.61 4.18c-2.81 2.79-4.21 6.75-4.21 11.87v.88c0 5.19 1.4 9.18 4.21 11.96 2.81 2.79 6.35 4.18 10.61 4.18M190.43 39.01c0-4.86-.97-9.11-2.91-12.74s-4.64-6.47-8.09-8.51c-3.46-2.04-7.48-3.07-12.06-3.07s-8.85 1.04-12.4 3.11a22.3 22.3 0 0 0-4.76 3.72c-.18.18-.35.36-.52.55-.51.56-.98 1.15-1.44 1.77-.6.82-1.16 1.7-1.67 2.62-2.03 3.7-3.05 8.04-3.05 13.04v1.17c0 5 1.02 9.34 3.05 13.04s4.86 6.57 8.48 8.61c3.61 2.04 7.85 3.07 12.69 3.07 4.26 0 7.82-.68 10.66-2.04s5.1-3.05 6.78-5.06q.705-.84 1.32-1.65c1.07-1.41 1.96-2.77 2.65-4.09l-8.33-4.28c-1.03 2.21-2.51 4.15-4.41 5.84-1.91 1.69-4.73 2.53-8.48 2.53-4.01 0-7.35-1.25-10.03-3.74-.17-.16-.33-.31-.48-.48-.16-.17-.31-.33-.46-.5-.36-.41-.69-.84-.99-1.3-1.3-1.96-2.05-4.3-2.24-7.01-.02-.26-.03-.52-.04-.78h36.72v-3.8Zm-30.04-3.89-5.1 5.12h-1.68c-.03-1.48 0-3.37.19-5.12.1-.67.22-1.31.39-1.93.16-.62.35-1.21.58-1.77.75-1.88 1.87-3.47 3.34-4.77 2.36-2.08 5.41-3.11 9.16-3.11s6.78 1.04 9.11 3.11c2.33 2.08 3.65 4.9 3.97 8.46h-19.96ZM303.6 39.01c0-4.86-.97-9.11-2.91-12.74s-4.64-6.47-8.09-8.51c-3.46-2.04-7.48-3.07-12.06-3.07s-8.85 1.04-12.4 3.11a22.3 22.3 0 0 0-4.76 3.72c-.17.18-.35.36-.52.55q-.765.84-1.44 1.77c-.6.82-1.16 1.7-1.67 2.62-2.03 3.7-3.05 8.04-3.05 13.04v1.17c0 5 1.02 9.34 3.05 13.04s4.86 6.57 8.48 8.61c3.61 2.04 7.85 3.07 12.69 3.07 4.26 0 7.81-.68 10.66-2.04 2.84-1.36 5.1-3.05 6.78-5.06q.705-.84 1.32-1.65c1.07-1.41 1.96-2.77 2.66-4.09l-8.33-4.28c-1.03 2.21-2.51 4.15-4.41 5.84-1.91 1.69-4.73 2.53-8.48 2.53-4.01 0-7.35-1.25-10.03-3.74-.17-.16-.33-.31-.48-.48-.16-.17-.31-.33-.46-.5-.36-.41-.69-.84-.99-1.3-1.3-1.96-2.05-4.3-2.24-7.01-.02-.26-.03-.52-.04-.78h36.73v-3.8Zm-30.04-3.89-5.1 5.12h-1.68c-.03-1.48 0-3.37.19-5.12.1-.67.22-1.31.39-1.93.16-.62.35-1.21.58-1.77.75-1.88 1.87-3.47 3.34-4.77 2.36-2.08 5.41-3.11 9.16-3.11s6.78 1.04 9.11 3.11c2.33 2.08 3.65 4.9 3.97 8.46h-19.96ZM312.21 65.37V14.69H322v5.64h1.55c.77-2.01 2.02-3.48 3.73-4.43s3.83-1.41 6.35-1.41h5.72v9.05h-6.1c-3.23 0-5.88.89-7.95 2.68s-3.1 4.52-3.1 8.22v30.93h-9.98ZM386.21 32.86c0-5.86-1.78-10.39-5.35-13.58s-8.47-4.79-14.7-4.79c-4.09 0-7.54.65-10.36 1.96s-5.09 3.03-6.81 5.18-2.97 4.5-3.75 7.03l9.34 3.03c.58-2.6 1.78-4.72 3.6-6.35s4.44-2.44 7.88-2.44 6.1.86 7.79 2.59c1.69 1.72 2.53 3.96 2.53 6.69v8.33h-1.62l-5.19-5.21h-7.79c-5.26 0-9.62 1.25-13.09 3.76s-5.21 6.2-5.21 11.09c0 3.26.8 6.02 2.38 8.3 1.59 2.28 3.73 4.01 6.42 5.18 2.7 1.17 5.76 1.76 9.19 1.76s5.97-.47 7.98-1.42c2.01-.94 3.54-2.05 4.57-3.32 1.04-1.27 1.81-2.43 2.34-3.47h1.56v8.21h8.27V32.85Zm-9.83 11.24c0 4.04-1.25 7.22-3.75 9.53-2.49 2.31-5.72 3.47-9.68 3.47-2.92 0-5.22-.65-6.91-1.96s-2.53-3.06-2.53-5.28.81-3.89 2.43-5.03c1.63-1.14 3.76-1.71 6.42-1.71h14.01v.98Z" style="fill:#160f26;stroke-width:0"></path></g></svg>',
    );
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
    $("#ai_provider_logo").html("");
    $("#mqc_anonymize_samples_switch").show();
    $("#ai_context_window_group").hide();
  } else if (providerId === "none") {
    $(".ai-generate-button-wrapper").hide();
    $(".ai-copy-button-wrapper").hide();
    $("#ai_provider_info").html("Remove 'Summarize' buttons from report if you're not going to use AI summaries.");
    $("#ai_provider_logo").html("");
    $("#ai_model_group").hide();
    $("#ai_api_key_group").hide();
    $("#ai_endpoint_group").hide();
    $("#ai_query_options_group").hide();
    $("#mqc_anonymize_samples_switch").hide();
    $("#ai_context_window_group").hide();
  } else {
    // Custom or other providers
    $(".ai-generate-button-wrapper").show();
    $(".ai-copy-button-wrapper").hide();
    $("#ai_model_group").show();
    $("#ai_api_key_group").show();
    $("#mqc_anonymize_samples_switch").show();

    if (providerId === "custom") {
      $("#ai_endpoint_group").show();
      $("#ai_query_options_group").show();
      $("#ai_context_window_group").show();
    } else {
      $("#ai_endpoint_group").hide();
      $("#ai_query_options_group").hide();
      $("#ai_context_window_group").hide();
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
function storeContextWindow(value) {
  saveToLocalStorage(`mqc_ai_context_window`, value.toString());
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
  $(".mqc_switch_wrapper.mqc_anonymize_samples").click(function (e) {
    e.preventDefault();
    const switchElem = $(this).find(".mqc_switch");
    const isEnabled = switchElem.hasClass("on");
    if (isEnabled) {
      switchElem.removeClass("on").addClass("off").text("off");
      storeSampleAnonymizationEnabled(false);
    } else {
      switchElem.removeClass("off").addClass("on").text("on");
      storeSampleAnonymizationEnabled(true);
    }
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
