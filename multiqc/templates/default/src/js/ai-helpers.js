import DOMPurify from "dompurify";

function handleStreamError(error) {
  if (error.toString().includes("Failed to fetch")) {
    error = "Failed to connect to AI provider. Please check your internet connection and the API key, and try again.";
  }
  return error;
}

function isReasoningModel(model) {
  if (!model) return false;
  const reasoningModels = [
    // OpenAI reasoning models
    "o1",
    "o1-preview",
    "o1-mini",
    "o3",
    "o3-mini",
    "o3-pro",
    "o4-mini",
    // Anthropic Claude 4 series (extended thinking models)
    "claude-3-7-sonnet-latest",
    "claude-sonnet-4-5",
    "claude-sonnet-4-0",
    "claude-haiku-4-5",
    "claude-haiku-4-0",
    "claude-opus-4-5",
    "claude-opus-4-0",
  ];
  return reasoningModels.some((prefix) => model.toLowerCase().startsWith(prefix));
}

function getStoredReasoningEffort() {
  return localStorage.getItem("ai_reasoning_effort") || "medium";
}

function getStoredMaxCompletionTokens() {
  const stored = localStorage.getItem("ai_max_completion_tokens");
  if (stored) {
    try {
      return parseInt(stored);
    } catch (e) {
      console.error("Error parsing stored max completion tokens", e);
    }
  }
  return 4000; // Default value
}

function getStoredThinkingBudgetTokens() {
  const stored = localStorage.getItem("ai_thinking_budget_tokens");
  if (stored) {
    try {
      return parseInt(stored);
    } catch (e) {
      console.error("Error parsing stored thinking budget tokens", e);
    }
  }
  return 10000; // Default value for Anthropic extended thinking
}

function getStoredExtendedThinking() {
  const stored = localStorage.getItem("ai_extended_thinking");
  return stored === "true"; // Default to false
}

function storeReasoningEffort(effort) {
  localStorage.setItem("ai_reasoning_effort", effort);
}

function storeMaxCompletionTokens(tokens) {
  localStorage.setItem("ai_max_completion_tokens", tokens.toString());
}

function storeThinkingBudgetTokens(tokens) {
  localStorage.setItem("ai_thinking_budget_tokens", tokens.toString());
}

function storeExtendedThinking(enabled) {
  localStorage.setItem("ai_extended_thinking", enabled.toString());
}

// Load reasoning model configuration on page load
$(function () {
  // Set up event handlers for reasoning model controls if they exist
  $("#ai-reasoning-effort").change(function () {
    storeReasoningEffort($(this).val());
  });

  $("#ai-max-completion-tokens").change(function () {
    storeMaxCompletionTokens($(this).val());
  });

  $("#ai-thinking-budget-tokens").change(function () {
    storeThinkingBudgetTokens($(this).val());
  });

  $("#ai-extended-thinking").change(function () {
    storeExtendedThinking($(this).is(":checked"));
  });

  // Load stored values into controls if they exist
  const storedEffort = getStoredReasoningEffort();
  if ($("#ai-reasoning-effort").length && storedEffort) {
    $("#ai-reasoning-effort").val(storedEffort);
  }

  const storedTokens = getStoredMaxCompletionTokens();
  if ($("#ai-max-completion-tokens").length && storedTokens) {
    $("#ai-max-completion-tokens").val(storedTokens);
  }

  const storedThinkingTokens = getStoredThinkingBudgetTokens();
  if ($("#ai-thinking-budget-tokens").length && storedThinkingTokens) {
    $("#ai-thinking-budget-tokens").val(storedThinkingTokens);
  }

  const storedExtendedThinking = getStoredExtendedThinking();
  if ($("#ai-extended-thinking").length) {
    $("#ai-extended-thinking").prop("checked", storedExtendedThinking);
  }

  // Show/hide reasoning model controls based on selected model
  $("#ai-model").change(function () {
    const modelName = $(this).val();
    const isReasoning = isReasoningModel(modelName);
    const isClaudeReasoning =
      isReasoning &&
      (modelName.toLowerCase().startsWith("claude-opus-4") ||
        modelName.toLowerCase().startsWith("claude-sonnet-4") ||
        modelName.toLowerCase().startsWith("claude-haiku-4"));
    const isOpenAIReasoning = isReasoning && !isClaudeReasoning;

    // Show/hide reasoning-specific controls if they exist
    $("#ai_reasoning_effort_group").toggle(isOpenAIReasoning);
    $("#ai_max_completion_tokens_group").toggle(isOpenAIReasoning);
    $("#ai_extended_thinking_group").toggle(isClaudeReasoning);
    $("#ai_thinking_budget_tokens_group").toggle(isClaudeReasoning);

    if (isReasoning) {
      console.log(
        `Reasoning model selected: ${modelName} (${
          isClaudeReasoning ? "Claude extended thinking" : "OpenAI reasoning"
        })`,
      );
    }
  });

  // Initial state based on current model
  const currentModel = $("#ai-model").val();
  if (currentModel) {
    const isReasoning = isReasoningModel(currentModel);
    const isClaudeReasoning =
      isReasoning &&
      (currentModel.toLowerCase().startsWith("claude-opus-4") ||
        currentModel.toLowerCase().startsWith("claude-sonnet-4") ||
        currentModel.toLowerCase().startsWith("claude-haiku-4"));
    const isOpenAIReasoning = isReasoning && !isClaudeReasoning;

    $("#ai_reasoning_effort_group").toggle(isOpenAIReasoning);
    $("#ai_max_completion_tokens_group").toggle(isOpenAIReasoning);
    $("#ai_extended_thinking_group").toggle(isClaudeReasoning);
    $("#ai_thinking_budget_tokens_group").toggle(isClaudeReasoning);
  }
});

window.runStreamGeneration = function ({
  onStreamStart,
  onStreamNewToken,
  onStreamError,
  onStreamComplete,
  systemPrompt,
  userMessage,
  tags = [],
  title = "",
}) {
  const providerId = $("#ai-provider").val();
  const provider = window.AI_PROVIDERS[providerId];
  const modelName = $("#ai-model").val();
  const apiKey = $("#ai-api-key").val();
  const customEndpoint = $("#ai-endpoint").val();

  const fetchOptions = {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
  };

  if (provider.name === "Seqera AI") {
    if (apiKey) {
      fetchOptions.headers.Authorization = `Bearer ${apiKey}`;
    }
    fetchOptions.body = JSON.stringify({
      message: title + "\n\n:::details\n\n" + systemPrompt + "\n\n" + userMessage + "\n\n:::\n\n",
      stream: true,
      tags: ["multiqc", ...tags],
      title: title,
    });

    fetch(`${seqeraApiUrl}/internal-ai/report-summary`, fetchOptions)
      .then((response) => {
        if (!response.ok) {
          return response.json().then((errorData) => {
            const error = `HTTP ${response.status}: ${response.statusText} ${
              errorData.error?.message || "Unknown error"
            }`;
            onStreamError(error);
            throw new Error(error);
          });
        }
        return response.body.getReader();
      })
      .then((reader) =>
        decodeStream(providerId, reader, onStreamStart, onStreamNewToken, onStreamError, onStreamComplete),
      )
      .catch((error) => {
        onStreamError(handleStreamError(error));
      });
  } else if (provider.name === "OpenAI" || provider.name === "Custom") {
    let body = {};
    if (provider.name === "Custom") {
      try {
        body = JSON.parse($("#ai-query-options").val());
      } catch (e) {
        console.error("Error parsing extra query options", e);
      }
    }

    // Check if this is a reasoning model
    const isReasoning = isReasoningModel(modelName);
    const isClaudeReasoning =
      isReasoning &&
      (modelName.toLowerCase().startsWith("claude-opus-4") ||
        modelName.toLowerCase().startsWith("claude-sonnet-4") ||
        modelName.toLowerCase().startsWith("claude-haiku-4"));

    if (isReasoning && !isClaudeReasoning) {
      console.log(`Using OpenAI reasoning model: ${modelName}`);

      // For OpenAI reasoning models, use developer messages and specific parameters
      const messages = [
        {
          role: "developer",
          content:
            "You are a helpful assistant expert in bioinformatics. Please provide clear, well-formatted responses using markdown where appropriate.",
        },
        { role: "user", content: systemPrompt + "\n\n" + userMessage },
      ];

      body = {
        ...body,
        model: modelName,
        messages: messages,
        stream: true,
        max_completion_tokens: getStoredMaxCompletionTokens(),
        reasoning_effort: getStoredReasoningEffort(),
      };

      // Remove unsupported parameters for reasoning models
      delete body.temperature;
      delete body.top_p;
      delete body.presence_penalty;
      delete body.frequency_penalty;
      delete body.max_tokens; // Use max_completion_tokens instead

      console.log(
        `OpenAI reasoning model parameters: max_completion_tokens=${body.max_completion_tokens}, reasoning_effort=${body.reasoning_effort}`,
      );
    } else if (isClaudeReasoning) {
      const extendedThinkingEnabled = getStoredExtendedThinking();

      if (extendedThinkingEnabled) {
        console.log(`Using Claude 4 model with extended thinking: ${modelName}`);

        const thinkingBudgetTokens = getStoredThinkingBudgetTokens();

        body = {
          ...body,
          model: modelName,
          messages: [{ role: "user", content: systemPrompt + "\n\n" + userMessage }],
          stream: true,
          thinking: {
            type: "enabled",
            budget_tokens: thinkingBudgetTokens,
          },
        };

        console.log(`Extended thinking parameters: budget_tokens=${thinkingBudgetTokens}`);
      } else {
        console.log(`Using Claude 4 model without extended thinking: ${modelName}`);

        body = {
          ...body,
          model: modelName,
          messages: [{ role: "user", content: systemPrompt + "\n\n" + userMessage }],
          stream: true,
        };
      }
    } else {
      // Regular models
      body = {
        ...body,
        model: modelName,
        messages: [{ role: "user", content: systemPrompt + "\n\n" + userMessage }],
        stream: true,
      };
    }

    fetchOptions.body = JSON.stringify(body);
    fetchOptions.headers.Authorization = `Bearer ${apiKey}`;

    const endpoint = provider.name === "Custom" ? customEndpoint : "https://api.openai.com/v1/chat/completions";
    fetch(endpoint, fetchOptions)
      .then(async (response) => {
        if (!response.ok) {
          const errorData = await response.json();
          onStreamError(
            `HTTP ${response.status}: ${response.statusText} ${errorData.error?.message || "Unknown error"}`,
          );
          throw new Error(errorData.error?.message || "Unknown error");
        }
        return response.body.getReader();
      })
      .then((reader) =>
        decodeStream(providerId, reader, onStreamStart, onStreamNewToken, onStreamError, onStreamComplete),
      )
      .catch((error) => {
        onStreamError(handleStreamError(error));
      });
  } else if (provider.name === "Anthropic") {
    fetchOptions.headers = {
      ...fetchOptions.headers,
      "x-api-key": apiKey,
      "anthropic-version": "2023-06-01",
      "anthropic-dangerous-direct-browser-access": "true",
    };
    fetchOptions.body = JSON.stringify({
      model: modelName,
      max_tokens: 4096,
      messages: [{ role: "user", content: systemPrompt + "\n\n" + userMessage }],
      stream: true,
    });

    fetch("https://api.anthropic.com/v1/messages", fetchOptions)
      .then(async (response) => {
        if (!response.ok) {
          const errorData = await response.json();
          const errMessage = errorData.error ? `${errorData.error.type}: ${errorData.error.message}` : "Unknown error";
          onStreamError(`HTTP ${response.status}: ${response.statusText} ${errMessage}`);
          throw new Error(errMessage);
        }
        return response.body.getReader();
      })
      .then((reader) =>
        decodeStream(providerId, reader, onStreamStart, onStreamNewToken, onStreamError, onStreamComplete),
      )
      .catch((error) => {
        onStreamError(handleStreamError(error));
      });
  } else {
    onStreamError(`Unsupported AI provider: ${provider.name}`);
  }
};

function decodeStream(providerId, reader, onStreamStart, onStreamNewToken, onStreamError, onStreamComplete) {
  // Decode the stream. Supports Anthropic and OpenAI streaming responses
  const decoder = new TextDecoder();
  let buffer = "";

  let model = undefined;
  let threadId = undefined;
  let started = false;

  return recursivelyProcessStream();
  function recursivelyProcessStream() {
    // Inner function to recursively process the stream reader
    return reader
      .read()
      .then(({ value, done }) => {
        if (done) {
          onStreamComplete(threadId);
          return;
        }

        buffer += decoder.decode(value, { stream: true });
        const lines = buffer.split("\n");

        // Process all complete lines
        buffer = lines.reduce((acc, line) => {
          line = line.trim();
          if (!line) return acc;
          if (!line.startsWith("data: ")) {
            if (line.includes('"type":"invalid_request_error"')) {
              // "{"type":"error","error":{"type":"invalid_request_error","message":"messages.1.content: Input should be a valid list"}}"
              onStreamError(line);
              return acc;
            }
            return acc;
          }
          let jsonString = line.slice(6);
          if (jsonString === "[DONE]") {
            // OpenAI last line is "data: [DONE]"
            onStreamComplete();
            return acc;
          }
          let data;
          try {
            data = JSON.parse(jsonString);
          } catch (e) {
            // Unexpected JSON format. Ignore the line, but if somethings is really wrong,
            // we should anyway call onStreamComplete in the end
            // onStreamError(`Error parsing JSON from streaming response. JSON: ${jsonString}, error: ${e}`);
            return acc;
          }

          let type = undefined;
          let content = undefined;
          let role = undefined;
          let finishReason = undefined;

          // Detect provider
          if (providerId === "seqera") {
            type = data.type; // on_chat_model_start, on_chat_model_stream, on_chat_model_end
            content = data.content;
            model = data.metadata?.ls_model_name;
            threadId = threadId || data.thread_id;
            finishReason = data.response_metadata?.stop_reason;
            if (finishReason && finishReason !== "end_turn") {
              type = "error";
              error = `Streaming unexpectedly stopped. Reason: ${finishReason}`;
            }
          } else if (providerId === "openai" || providerId === "custom") {
            content = data.choices[0].delta.content;
            model = data.model;
            role = data.choices[0].delta.role;
            finishReason = data.choices[0].finish_reason;
            // OpenAI doesn't define type, so we need to infer it from the other fields
            if (content) {
              type = "on_chat_model_stream";
            } else if (finishReason === "stop") {
              type = "on_chat_model_end";
            } else if (role === "assistant") {
              type = "on_chat_model_start";
            } else if (finishReason && finishReason !== "stop") {
              type = "error";
              error = `Streaming unexpectedly stopped. Reason: ${finishReason}`;
            } else {
              type = "unknown";
            }
          } else if (providerId === "anthropic") {
            type = data.type;
            if (type === "content_block_delta" && data.delta.type === "text_delta") {
              content = data.delta.text;
              type = "on_chat_model_stream";
            } else if (type === "message_start") {
              model = data.message.model;
              role = data.message.role;
              type = "on_chat_model_start";
            } else if (data.finish_reason) {
              if (data.finish_reason === "end_turn") {
                type = "on_chat_model_end";
              } else {
                type = "error";
                error = `Streaming unexpectedly stopped. Reason: ${data.finish_reason}`;
                if (data.error) error += `. Error: ${data.error}`;
              }
            } else if (data.error) {
              type = "error";
              error = `Streaming unexpectedly stopped. Error: ${data.error}`;
            }
          }

          // Handle different event types
          switch (type) {
            case "on_chat_model_start":
              started = true;
              onStreamStart(model);
              break;
            case "on_chat_model_stream":
              if (!started) {
                started = true;
                onStreamStart(model);
              }
              if (content) onStreamNewToken(content);
              break;
            case "on_chat_model_end":
              onStreamComplete(threadId);
              return acc;
            case "error":
              onStreamError(error);
              break;
          }
          return acc;
        }, "");

        return recursivelyProcessStream();
      })
      .catch((error) => onStreamError(error));
  }
}

function deanonymizeSampleNames(text) {
  // Convert pseudonyms back to original sample names in the text.
  // Only applies when anonymization is enabled and pseudonym map exists
  if (!aiPseudonymMap) return text;
  if (!getStoredSampleAnonymizationEnabled()) return text;

  // Create reverse mapping from pseudonym to original name
  const reverseMap = Object.fromEntries(
    Object.entries(aiPseudonymMap).map(([original, pseudonym]) => [pseudonym, original]),
  );

  // Replace pseudonyms with original names, being careful to only replace whole words
  // and preserve the directive syntax
  for (const [pseudonym, original] of Object.entries(reverseMap)) {
    // Look for pseudonym in :sample[SAMPLE_1]{.text-*} directives
    text = text.replace(
      new RegExp(`:sample\\[${escapeRegExp(pseudonym)}\\](\\{[^}]+\\})`, "g"),
      `:sample[${original}]$1`,
    );

    // Look for standalone pseudonyms (not in directives)
    text = text.replace(new RegExp(`\\b${escapeRegExp(pseudonym)}\\b`, "g"), original);
  }

  return text;
}

// Helper function to escape special characters in string for use in RegExp
function escapeRegExp(string) {
  return string.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

window.markdownToHtml = function (text) {
  if (!text) return "";

  text = deanonymizeSampleNames(text);

  // Markdown first, so that code spans escape their own angle brackets and "`x < 5`"
  // does not come out double-escaped.
  try {
    let converter = new showdown.Converter({
      literalMidWordUnderscores: true, // Prevents interpretation of underscores within words
    });
    text = converter.makeHtml(text);
  } catch (e) {
    // Fall through and sanitise the unconverted text
  }

  // The text can carry sample names straight from the analysed files, so a crafted name
  // could otherwise inject markup into the report. Default profile: ordinary formatting
  // and links survive, scripts, event handlers and javascript: URLs do not.
  // The markdown -> sanitise -> directives order matters, and this function is forked
  // into multiqc/templates/original/assets/js/ai-helpers.js - keep the two in step.
  text = DOMPurify.sanitize(text);

  // Directives last: this is our own markup, and the captured text is already
  // sanitised, so re-inserting it here is safe. Do not move this above the sanitise.
  // :span[text]{.text-color} -> <span class="text-color">... (preserving underscores)
  text = text.replace(/:span\[([^\]]+?)\]\{\.text-(green|red|yellow)\}/g, (match, p1, p2) => {
    return `<span class="text-${p2}">${p1}</span>`;
  });

  // :sample[text]{.text-color} -> <sample... (preserving underscores)
  text = text.replace(/:sample\[([^\]]+?)\]\{\.text-(green|red|yellow)\}/g, (match, p1, p2) => {
    return `<sample data-bs-toggle="tooltip" title="Click to highlight in the report" class="text-${p2}">${p1}</sample>`;
  });

  return text;
};

const missingPromptFallback = "Analyze the supplied MultiQC data and summarize important quality-control findings.";
const reportPrompts = typeof aiPrompts === "undefined" ? {} : aiPrompts;

window.systemPromptReportShort = reportPrompts.short || missingPromptFallback;
window.systemPromptReportFull = reportPrompts.full || missingPromptFallback;
window.systemPromptPlot = reportPrompts.section || missingPromptFallback;
window.aiSummaryFull = Boolean(reportPrompts.full_summary);
