async function streamGeneration(
  onStreamStart,
  onStreamNewToken,
  onStreamError,
  onStreamComplete,
  systemPrompt,
  userMessage,
  tags = [],
) {
  function handleAnthropicOrOpenAiStream(reader) {
    const decoder = new TextDecoder();
    let buffer = "";

    function processStream() {
      return reader
        .read()
        .then(({ value, done }) => {
          if (done) {
            onStreamComplete();
            return;
          }

          buffer += decoder.decode(value, { stream: true });
          const lines = buffer.split("\n");

          // Process all complete lines
          buffer = lines.reduce((acc, line) => {
            line = line.trim();
            if (!line) return acc;
            if (line.startsWith("data: ")) {
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
                console.log(`Error parsing JSON from streaming response. JSON: ${jsonString}, error: ${e}`);
                return acc;
              }
              let type = data.type;
              let content = undefined;
              let model = undefined;
              let role = undefined;
              let finish_reason = undefined;

              // OpenAI doesn't define type
              if (type === undefined) {
                content = data.choices[0].delta.content;
                model = data.model;
                role = data.choices[0].delta.role;
                finish_reason = data.choices[0].finish_reason;
                if (role === "assistant") {
                  type = "message_start";
                } else if (finish_reason === "stop") {
                  type = "message_stop";
                } else if (content) {
                  type = "content_block_delta";
                } else if (finish_reason && finish_reason !== "stop") {
                  type = "error";
                } else {
                  type = "unknown";
                }
                error = finish_reason;
              } else {
                if (type === "content_block_delta" && data.delta.type === "text_delta") content = data.delta.text;
                if (type === "message_start") {
                  model = data.message.model;
                  role = data.message.role;
                }
                finish_reason = data.finish_reason;
                error = data.error;
              }

              console.log(type, line);

              // Handle different event types
              switch (type) {
                case "message_start":
                  onStreamStart(model);
                  break;
                case "content_block_delta":
                  if (content) onStreamNewToken(content);
                  break;
                case "message_stop":
                  onStreamComplete();
                  return acc;
                case "error":
                  onStreamError(error);
                  break;
              }
            }
            return acc;
          }, "");

          return processStream();
        })
        .catch((error) => onStreamError(error));
    }

    try {
      return processStream();
    } catch (e) {
      onStreamError(e);
    }
  }

  const provider = getStoredProvider();
  const modelName = getStoredModelName(provider);
  const apiKey = getStoredApiKey(provider);

  if (provider === "Seqera AI") {
    response = await fetch(`${seqeraApiUrl}/invoke`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        ...(apiKey ? { Authorization: `Bearer ${apiKey}` } : {}),
      },
      body: JSON.stringify({
        system_message: systemPrompt,
        user_message: userMessage,
        model: modelName,
        tags: ["multiqc", ...tags],
        stream: true,
        cache: false,
      }),
    });
  } else if (provider === "OpenAI") {
    fetch(`https://api.openai.com/v1/chat/completions`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        Authorization: `Bearer ${apiKey}`,
      },
      body: JSON.stringify({
        model: modelName,
        messages: [
          { role: "user", content: systemPrompt },
          { role: "user", content: userMessage },
        ],
        stream: true,
      }),
    })
      .then((response) => response.body.getReader())
      .then(handleAnthropicOrOpenAiStream)
      .catch((error) => onStreamError(error));
  } else if (provider === "Anthropic") {
    fetch(`https://api.anthropic.com/v1/messages`, {
      method: "POST",
      headers: {
        "content-type": "application/json",
        "x-api-key": apiKey,
        "anthropic-version": "2023-06-01",
        "anthropic-dangerous-direct-browser-access": "true",
      },
      body: JSON.stringify({
        model: modelName,
        max_tokens: 4096,
        messages: [
          { role: "user", content: systemPrompt },
          { role: "user", content: userMessage },
        ],
        stream: true,
      }),
    })
      .then((response) => response.body.getReader())
      .then(handleAnthropicOrOpenAiStream)
      .catch((error) => onStreamError(error));
  } else {
    onStreamError(`Unsupported AI provider: ${provider}`);
  }
}

function markdownToHtml(text) {
  if (!text) return "";

  // Convert directives :span[text]{.text-color} -> <span>... (preserving underscores)
  text = text.replace(/:span\[([^\]]+?)\]\{\.text-(green|red|yellow)\}/g, (match, p1, p2) => {
    return `<span class="text-${p2}">${p1}</span>`;
  });

  // Convert directives :sample[text]{.text-color} -> <sample>... (preserving underscores)
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

let systemPromptReport = `\
You are an expert in bioinformatics, sequencing technologies, genomics data analysis, and adjacent fields.

You are given findings from a MultiQC report, generated by a bioinformatics workflow.
MultiQC consists of so called modules that support different tools that output QC metrics
(e.g. bclconvert, FastQC, samtools stats, bcftools stats, fastp, Picard, SnpEff, etc), and
it aggregates results from different tools. It outputs a "General Statistics" section that
has a table with a summary of key metrics from all modules across each sample. That table
is followed by module-specific sections that usually have more detail on the same samples.

You are given data from such a report, split by section. Your task is to analyse this data and
generate an overall summary for the results. Please don't print any introductory words, just get to the 
point. You task is to just generate a concise summary of the report, nothing else. 
Don't waste words: mention only the important QC issues. If there are no issues, just say so. 
Try to format the response with bullet points. Follow up with recommendations for the next steps.

Make sure to use markdown to format your reponse for readability. Use directives with pre-defined classes
.text-green, .text-red, and .text-yellow to highlight severity, e.g. :span[39.2%]{.text-red}. 
If there are any sample names mentioned, or sample name prefixes or suffixes, you must warp them in
a sample directive, making sure to use same color classes as for severity, for example: :sample[A1001.2003]{.text-yellow}
or :sample[A1001]{.text-yellow}. But never put multiple sample names inside one directive.

Please do not add any extra headers to the response.

You must use only multiples of 4 spaces to indent nested lists.

This is the example response:

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

let systemPromptPlot = `\
You are an expert in bioinformatics, sequencing technologies, genomics data analysis, and adjacent fields.

You are given plot data from a MultiQC report that was generated by a bioinformatics workflow.
Your task is to analyse the data, and give a very short and concise overall summary of the results.
Don't waste words: mention only the important QC issues. If there are no issues, just say so.
Limit it to 1-2 sentences.

Make sure to use markdown to format your reponse for readability. Use directives with pre-defined classes
.text-green, .text-red, and .text-yellow to highlight severity, e.g. :span[39.2%]{.text-red}. 
If there are any sample names mentioned, or sample name prefixes or suffixes, you must warp them in
a sample directive, making sure to use same color classes as for severity, for example: :sample[A1001.2003]{.text-yellow}
or :sample[A1001]{.text-yellow}. But never put multiple sample names inside one directive.

Please do not add any extra headers to the response.

Make sure to use a multiple of 4 spaces to indent nested lists.`;
