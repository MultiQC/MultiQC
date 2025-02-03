function handleStreamError(error) {
  if (error.toString().includes("Failed to fetch")) {
    error = "Failed to connect to AI provider. Please check your internet connection and the API key, and try again.";
  }
  return error;
}

function runStreamGeneration({
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
  const provider = AI_PROVIDERS[providerId];
  const modelName = $("#ai-model").val();
  const apiKey = $("#ai-api-key").val();

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

    fetch(`${seqeraApiUrl}/internal-ai/query`, fetchOptions)
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
  } else if (provider.name === "OpenAI") {
    fetchOptions.headers.Authorization = `Bearer ${apiKey}`;
    fetchOptions.body = JSON.stringify({
      model: modelName,
      messages: [
        { role: "user", content: systemPrompt },
        { role: "user", content: userMessage },
      ],
      stream: true,
    });

    fetch("https://api.openai.com/v1/chat/completions", fetchOptions)
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
      messages: [
        { role: "user", content: systemPrompt },
        { role: "user", content: userMessage },
      ],
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
}

function decodeStream(providerId, reader, onStreamStart, onStreamNewToken, onStreamError, onStreamComplete) {
  // Decode the stream. Supports Anthropic and OpenAI streaming responses
  const decoder = new TextDecoder();
  let buffer = "";

  let model = undefined;
  let threadId = undefined;

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
            onStreamError(`Error parsing JSON from streaming response. JSON: ${jsonString}, error: ${e}`);
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
          } else if (providerId === "openai") {
            content = data.choices[0].delta.content;
            model = data.model;
            role = data.choices[0].delta.role;
            finishReason = data.choices[0].finish_reason;
            // OpenAI doesn't define type, so we need to infer it from the other fields
            if (role === "assistant") {
              type = "on_chat_model_start";
            } else if (finishReason === "stop") {
              type = "on_chat_model_end";
            } else if (content) {
              type = "on_chat_model_stream";
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
              onStreamStart(model);
              break;
            case "on_chat_model_stream":
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

function markdownToHtml(text) {
  if (!text) return "";

  text = deanonymizeSampleNames(text);

  // Convert directives :span[text]{.text-color} -> <span class="text-color">... (preserving underscores)
  text = text.replace(/:span\[([^\]]+?)\]\{\.text-(green|red|yellow)\}/g, (match, p1, p2) => {
    return `<span class="text-${p2}">${p1}</span>`;
  });

  // Convert directives :sample[text]{.text-color} -> <sample... (preserving underscores)
  text = text.replace(/:sample\[([^\]]+?)\]\{\.text-(green|red|yellow)\}/g, (match, p1, p2) => {
    return `<sample data-toggle="tooltip" title="Click to highlight in the report" class="text-${p2}">${p1}</sample>`;
  });

  // Convert markdown to html
  try {
    let converter = new showdown.Converter({
      literalMidWordUnderscores: true, // Prevents interpretation of underscores within words
    });
    return converter.makeHtml(text);
  } catch (e) {
    return text;
  }
}

let multiqcDescription = `\
You are an expert in bioinformatics, sequencing technologies, genomics data analysis, and adjacent fields.

You are given findings from a MultiQC report, generated by a bioinformatics workflow.
MultiQC supports various bioinformatics tools that output QC metrics, and aggregates those metrics
into a single report. It outputs a "General Statistics" table with key metrics for each sample across
all tools. That table is followed by more detailed sections from specific tools, that can include tables,
as well as plots of different types (bar plot, line plot, scatter plot, heatmap, etc.)
`;

let systemPromptReport =
  multiqcDescription +
  `
You are given data from such a report. Your task is to analyse this data and
generate an overall summary for the results.

Please don't print any introductory words, just get to the point.
You task is to just generate a concise summary of the report, nothing else.
Don't waste words: mention only the important QC issues. If there are no issues, just say so.
Try to format the response with bullet points. Please do not add any extra headers to the response.

Use markdown to format your reponse for readability. Use directives with pre-defined classes
.text-green, .text-red, and .text-yellow to highlight severity, e.g. :span[39.2%]{.text-red}.
If there are any sample names mentioned, or sample name prefixes or suffixes, you must warp them in
a sample directive, making sure to use same color classes as for severity, for example: :sample[A1001.2003]{.text-yellow}
or :sample[A1001]{.text-yellow}. But never put multiple sample names inside one directive.

You must use only multiples of 4 spaces to indent nested lists.
`;

let systemPromptReportShort =
  systemPromptReport +
  `
Limit the response to 1-2 bullet points. Two such examples of short summaries:

- :span[11/13 samples]{.text-green} show consistent metrics within expected ranges.
- :sample[A1001.2003]{.text-red} and :sample[A1001.2004]{.text-red} exhibit extremely high percentage of :span[duplicates]{.text-red} (:span[65.54%]{.text-red} and :span[83.14%]{.text-red}, respectively).

- All samples show good quality metrics with :span[75.7-77.0%]{.text-green} CpG methylation and :span[76.3-86.0%]{.text-green} alignment rates
- :sample[2wk]{.text-yellow} samples show slightly higher duplication (:span[11-15%]{.text-yellow}) compared to :sample[1wk]{.text-green} samples (:span[6-9%]{.text-green})'
`;

let systemPromptReportFull =
  systemPromptReport +
  `
Follow up with recommendations for the next steps.

This is the example response:

**Analysis**

- :sample[A1002]{.text-yellow} and :sample[A1003]{.text-yellow} groups (:span[11/13 samples]{.text-green}) show good quality metrics, with consistent GC content (38-39%), read lengths (125 bp), and acceptable levels of duplicates and valid pairs.
- :sample[A1001.2003]{.text-red} and :sample[A1001.2004]{.text-red} show severe quality issues:
    - Extremely high duplicate rates (:span[65.54%]{.text-red} and :span[83.14%]{.text-red})
    - Low percentages of valid pairs (:span[37.2%]{.text-red} and :span[39.2%]{.text-red})
    - High percentages of failed modules in FastQC (:span[33.33%]{.text-red})
    - Significantly higher total sequence counts (:span[141.9M]{.text-red} and :span[178.0M]{.text-red}) compared to other samples
    - FastQC results indicate that :sample[A1001.2003]{.text-red} and :sample[A1001.2004]{.text-red} have a slight :span[GC content]{.text-red} bias at 39.5% against most other samples having 38.0%, which indicates a potential contamination that could be the source of other anomalies in quality metrics.

- :sample[A1002-1007]{.text-yellow} shows some quality concerns:
    - Low percentage of valid pairs (:span[48.08%]{.text-yellow})
    - Low percentage of passed Di-Tags (:span[22.51%]{.text-yellow})

- Overrepresented sequences analysis reveals adapter contamination in several samples, particularly in :sample[A1001.2003]{.text-yellow} (up to :span[35.82%]{.text-yellow} in Read 1).
- HiCUP analysis shows that most samples have acceptable levels of valid pairs, with :sample[A1003]{.text-green} group generally performing better than :sample[A1002]{.text-yellow} group.

**Recommendations**

- Remove :sample[A1001.2003]{.text-red} and :sample[A1200.2004]{.text-red} from further analysis due to severe quality issues.
- Investigate the cause of low valid pairs and passed Di-Tags in :sample[A1002-1007]{.text-yellow}. Consider removing it if the issue cannot be resolved.
- Perform adapter trimming on all samples, particularly focusing on :sample[A1001]{.text-red} group.
- Re-run the Hi-C analysis pipeline after removing problematic samples and performing adapter trimming.
- Investigate the cause of higher duplication rates in :sample[A1002]{.text-yellow} group compared to :sample[A1003]{.text-green} group, although they are still within acceptable ranges.
- Consider adjusting the Hi-C protocol or library preparation steps to improve the percentage of valid pairs, especially for :sample[A1002]{.text-yellow} group.
`;

let systemPromptPlot =
  multiqcDescription +
  `
You are given a single MultiQC report section with a plot or a table.
Your task is to analyse the data and give a very short and concise overall summary of the results.
Don't waste words: mention only the important QC issues. If there are no issues, just say so.
Limit it to 1-2 sentences.

Make sure to use markdown to format your reponse for readability. Use directives with pre-defined classes
.text-green, .text-red, and .text-yellow to highlight severity, e.g. :span[39.2%]{.text-red}.
If there are any sample names mentioned, or sample name prefixes or suffixes, you must warp them in
a sample directive, making sure to use same color classes as for severity, for example: :sample[A1001.2003]{.text-yellow}
or :sample[A1001]{.text-yellow}. But never put multiple sample names inside one directive.

Please do not add any extra headers to the response.

Make sure to use a multiple of 4 spaces to indent nested lists.`;
