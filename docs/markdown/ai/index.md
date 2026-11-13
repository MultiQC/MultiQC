---
title: AI summaries
description: Using AI to summarise MultiQC reports
---

# AI summaries

MultiQC v1.27 and newer can generate AI-powered summaries of your reports. These can be created at two points:

- When creating the report, baked into the report HTML.
- Dynamically in the browser, while viewing an existing HTML report.

The AI summaries are generated using LLMs (large-language models) AI, with the following supported providers:

- [Seqera Co-Scientist](https://ai.seqera.io/)
- [OpenAI](https://openai.com/)
- [Anthropic](https://www.anthropic.com/)
- [AWS Bedrock](https://aws.amazon.com/bedrock/)
- MultiQC reports also have an option to copy a prompt to your clipboard, to paste into any provider you have access to

:::warning

Never rely on AI-generated summaries. Whilst these summaries can be useful to get you started quickly with a report, they may give inaccurate analysis and miss important details.

:::

AI summaries work by sending report data to an LLM provider of your choice, via an API over the internet.
Be aware of what data you are sending, and to who.
For more information, see [Seqera AI: Your privacy matters](https://seqera.io/ai-trust/).

## Choosing a provider

To use native summary generation, MultiQC needs to communicate with an LLM provider's API.
All supported services require an API key or credentials to work.
Remember: Treat your API keys like passwords and do not share them.

- [Seqera Co-Scientist](https://ai.seqera.io/)
  - Register for free at [seqera.io](https://seqera.io/)
  - Create a new key on the **Access tokens** page: [https://cloud.seqera.io/tokens](https://cloud.seqera.io/tokens)
- [OpenAI](https://openai.com/)
  - Register at [platform.openai.com](https://platform.openai.com/signup) (NB: different to ChatGPT)
  - Add a payment method to your account to enable API usage beyond any trial credits
  - Create a new secret key on the _API Keys_ section [under your profile](https://platform.openai.com/api-keys)
- [Anthropic](https://www.anthropic.com/)
  - Sign up at [https://console.anthropic.com](https://console.anthropic.com)
  - Add a payment method to enable API access
  - Create a new key on on the _API Keys_ section in your [account settings](https://console.anthropic.com/settings/keys)
- [AWS Bedrock](https://aws.amazon.com/bedrock/)
  - Bedrock supports a plethora of models across many providers
  - Sign up, access credentials and payment are handled via an AWS account
- Other providers, via custom endpoint
  - Works with any OpenAI-compatible or Claude-compatible API, or your own Python class. See [Custom providers](#custom-providers) for details
- Other providers, via your clipboard
  - You can use buttons in MultiQC reports to copy a prompt to your clipboard,
    in order to manually summarise report data.
    See [Copying prompts](#copying-prompts) for instructions.

Seqera Co-Scientist is free to use.[^seqera-ai-usage-limits]
Use of other third-party APIs are billed by their respective providers based on consumption.
Seqera Co-Scientist uses the latest AI provider models under the hood.

### Choosing a model

If you're using OpenAI, Anthropic, AWS Bedrock or a custom endpoint you can choose the exact model used for report summaries.
This is done by setting `ai_model` in the MultiQC config.

- Anthropic: default `claude-opus-5-5`. Model names must contain `claude`.
  See the [Anthropic docs](https://platform.claude.com/docs/en/about-claude/models/overview).
- OpenAI: default `gpt-6.1-sol`.
  See the [OpenAI docs](https://developers.openai.com/api/docs/models).
- AWS Bedrock: required. Model names must be valid inputs to the `modelId` parameter of the `InvokeModel` API ([docs](https://docs.aws.amazon.com/bedrock/latest/APIReference/API_runtime_InvokeModel.html#API_runtime_InvokeModel_RequestSyntax)), for example `us.anthropic.claude-opus-5-5`.
- Custom endpoints: required. Use whatever name the endpoint expects.

This model is used during report generation and also set as the default toolbox panel setting for browser report summaries.

MultiQC doesn't keep a list of supported models. It reads the model family and version from the model name
(for example `claude-opus-5-5`, `gpt-6.1-sol`, or a Bedrock ID such as `us.anthropic.claude-sonnet-4-5-20250929-v1:0`)
and adapts the request to it, so newly released models work without a MultiQC update.
It also sends only what is needed: no sampling settings such as `temperature`, and no reasoning settings
unless you configure them, so each model runs with its own defaults.
If a model rejects a setting anyway, MultiQC sends the request again without it and prints a warning.

### Reasoning and thinking

Current models reason before they answer. You can tune this with the following options, which apply to all providers:

```yaml
# multiqc_config.yaml
ai_reasoning_effort: high # passed to the provider as given: OpenAI reasoning_effort, Claude output_config.effort
ai_max_completion_tokens: 32000 # maximum output tokens, including reasoning or thinking
```

- **`ai_reasoning_effort`**: the levels depend on the model, for example `low`, `medium`, `high`, `xhigh` or `max`.
  Lower effort is faster and cheaper. Not set by default, so each model uses its own default.
- **`ai_max_completion_tokens`**: defaults to 16,000 for Claude models. Other providers use the model's own maximum.
  Reasoning tokens count towards this limit, so setting it too low can leave no room for the answer.

For Claude models, you can also turn on extended thinking:

```yaml
ai_extended_thinking: true
ai_thinking_budget_tokens: 15000 # only used by Claude models older than 4.6
```

On Claude 4.6 and later this uses adaptive thinking, where the model decides how much to think.
Older models take a fixed token budget. When `ai_extended_thinking` is off, each model uses its own default:
the newest Claude models think adaptively by default, and `ai_reasoning_effort: low` keeps that short.

Reasoning and thinking tokens are billed by the provider, but are not shown in the summary.

### Refusals

Claude Opus, Sonnet and Fable models from version 5 run safety classifiers that can occasionally decline
legitimate requests, including some about biology. When using the Anthropic API directly,
MultiQC asks Anthropic to re-run a declined request on another model automatically.
If a request is still declined, or is declined on another provider, the report shows a short note instead of a summary.

## Summaries during report generation

MultiQC can generate AI summaries at run time, when generating reports.
Summary text is included within the report HTML as static text and will be visible to anyone viewing the report,
even when shared.

### MultiQC configuration

AI summaries are disabled by default when running MultiQC.
To generate them, you must enable them either on the command line or via a MultiQC config file.

- Command line flags:
  - `--ai` / `--ai-summary`: Generate a short report summary and put it on top of the report (fast)
  - `--ai-summary-full`: Generate a detailed version of the summary with analysis and recommendations (slower)
  - `--ai-provider <provider>`: Choose AI provider. One of `seqera`, `openai`, `anthropic`, `aws_bedrock` or `custom`. Default `seqera`
  - `--no-ai`: Disable AI toolbox and buttons in the report

- Alternatively, MultiQC configuration file:
  ```yaml
  ai_summary: false # Set to true for short summaries
  ai_summary_full: false # Set to true for  long summaries
  ai_provider: "seqera" # 'seqera', 'openai', 'anthropic', 'aws_bedrock' or 'custom'. Default: 'seqera'
  no_ai: false # Set to true to disable AI toolbox and buttons in the report
  ```

### Environment variables

You must also set your provider's API key in an environment variable in order to access its service
_(see [Choosing a provider](#choosing-a-provider) for how to get an API key)_.
API keys are supplied by setting the following environment variables in your shell:

```bash
export SEQERA_ACCESS_TOKEN="..."  # or TOWER_ACCESS_TOKEN
export OPENAI_API_KEY="..."
export ANTHROPIC_API_KEY="..."
```

It's possible to save these in an `.env` file instead of exporting to your shell's environment.
This `.env` file can either be in the current working directory or the MultiQC source code directory.
MultiQC uses the [python-dotenv](https://saurabh-kumar.com/python-dotenv/) package.

For AWS Bedrock, the client uses the [default `boto3` credential chain](https://boto3.amazonaws.com/v1/documentation/api/latest/guide/credentials.html).

If you run MultiQC without the appropriate key you will get a warning printed to the console,
but report generation will otherwise proceed without the summary. MultiQC will not return an error exit code.

:::tip

MultiQC configuration options can also be set using environment variables
(see [Config with environment variables](../getting_started/config.md#config-with-environment-variables)),
so you can set up everything, including the command line flags / config this way:

```bash
export MULTIQC_AI_SUMMARY=1
export SEQERA_ACCESS_TOKEN="..."
```

:::

Environment variables will only be used for `--ai-summary`/`--ai-summary-full` generation.
They are not saved by MultiQC and cannot be used for in-browser summary generation, within reports.

## In-browser AI summaries

In addition to summaries during report generation, MultiQC can also create summaries dynamically in reports.
This can be useful as the person viewing a report is often different than the person who generated it.
Summaries can be generated on demand when needed.

The AI toolbox and **Summarize** buttons are shown by default in all reports. To prevent this, run MultiQC with the `--no-ai` flag.
This can also be done on a per-user basis by selecting **Remove AI buttons** in the **AI Provider** dropdown in the AI toolbox.

Summaries generated in reports are _ephemeral_ and are not saved in the HTML.
If you generate a summary and share the report then others will not see it.
MultiQC tries to save the summary response within your browser's [local storage](https://www.w3schools.com/html/html5_webstorage.asp)
so that it shows the next time you open the same report, but this process is imperfect and might not always work.

### Configuring the AI provider

<div style={{float:"right", width:200}}>

![Configure AI providers within a MultiQC report](../../../docs/images/ai_toolbox_icon.png)

</div>

When you first try to generate a summary in the browser, you must supply the LLM provider's API key.
Open the AI settings by clicking the icon in the toolbox:

Then, choose an AI provider and enter the relevant API key
_(see [Choosing a provider](#choosing-a-provider) for how to get an API key)_.

:::info[Important]

API keys are stored _only_ in your browser's [local storage](https://www.w3schools.com/html/html5_webstorage.asp)
and are not shared if you send the HTML report to someone else.
They are used to send report data directly to your AI provider of choice.

:::

![Enter a provider API key in the report toolbox](../../../docs/images/ai_toolbox_keys.png)

### Summarising the report

Once your provider API key is configured, click **Summarize report** to generate an overview summary of the entire report.

![Button to summarize a MultiQC report](../../../docs/images/ai_summarize_button.png)

The summary text is interactive: click an <u>underlined</u> sample name to highlight that sample throughout the report:

![Click underlined sample names in the summary to highlight them in the report](../../../docs/images/ai_highlight_samples.gif)

### Section-level summaries

Besides a global report-level AI summary, you can generate a summary for each plot or table separately using buttons next to each section:

![Summarize with AI buttons in a report](../../../docs/images/ai_summarize_buttons.png)

### Copying prompts

If you have access to an LLM that is not directly supported by MultiQC, you can copy the exact prompt
that MultiQC uses to your clipboard. This can be pasted into whatever LLM that you have access to.

To do this, select **Copy prompts** as the LLM provider in the report toolbox AI tab.
The **Summarize** buttons will then change to **Copy prompt** buttons and instead of injecting
summaries into the report HTML, will copy the LLM prompt to your clipboard.

### Remove AI buttons

If you're suffering from AI-overload and don't want to see the AI summary features in your reports,
you can disable them by selecting **Remove AI buttons** in the toolbox as an AI provider.
This will remove all **Summarize** buttons from the report.

This is done at user level and will be stored in the browser's local storage and applied to all
MultiQC reports that you open.
You can also use `--no-ai` when generating reports, which removes this functionality from the HTML for all users.

## Using `llms-full.txt`

MultiQC always saves the full prompt and response to `multiqc_data/llms-full.txt` file,
regardless of whether the summary was generated during report generation or in the browser.

This file can be used to debug or further analyse the AI summary generation process.
It can be used to directly copy the prompt into your clipboard and use it with external services, e.g. ones with a larger context window.

## Continue chat

If using Seqera AI as a provider, you can click the **Chat with Seqera AI** button to open the Seqera AI
chat interface in a new tab in order to ask further questions.
This button is shown alongside the report-level summary after it's generated.

![Chat with Seqera AI button location](../../../docs/images/ai_chat_with_seqera_ai_button.png)

If you are logged in to [seqera.io](https://seqera.io) with the same user that generated the report summary,
the chat history with the report prompt and AI summary response will be loaded allowing you to
continue straight on with more in-depth questions.

![Seqera AI with MultiQC report history](../../../docs/images/ai_continue_chat.png)

## Context window

A context window refers to the amount of text (in tokens) that an AI model can consider at once
when processing input and generating responses. Current models take 200k to 1M tokens or more.
MultiQC sends as much of the report as fits, as more data gives better summaries.

MultiQC finds the context window of the model from the first of these that gives an answer:

1. The `ai_custom_context_window` config option (or the **Context Window** field in the report toolbox, for custom endpoints).
2. The provider: the Anthropic models API, or the models listing of a custom endpoint
   (vLLM, OpenRouter, LM Studio and LiteLLM report context windows there).
3. What MultiQC knows about the model family, for Claude and OpenAI models.
4. 128,000 tokens, for models it doesn't recognise.

It then fills the prompt in order of priority:

1. The list of tools used in the report.
2. The general statistics table (without hidden-by-default columns, if the full table doesn't fit).
3. Every report section that still fits, in report order. Sections that don't fit are listed by name, so the model knows about them.

If the provider says the prompt is still too long, MultiQC reads the real limit from the error,
rebuilds a shorter prompt and tries again.
If even the general statistics table doesn't fit, no summary is generated.

:::note

Some local servers, such as Ollama, silently cut prompts that are longer than their context window
instead of returning an error. Ollama's default context window (`num_ctx`) is much smaller than 128,000 tokens,
so set `ai_custom_context_window` to match your server's settings.

:::

Sending a whole large report costs more: providers bill per token, and some charge more for prompts over a size
threshold. Set `ai_custom_context_window` to cap how much is sent.

If you're unable to generate an AI summary, you can try the following:

- Hide additional columns in the general statistics table (see [Hiding Columns](../reports/customisation.md#hiding-columns)).
- Hide General statistics data in the browser, and request the AI summary dynamically:
  - Hide columns with the **Configure columns** button
  - Filter shown samples dynamically with the toolbox
- Use a model with a larger context window.
- Copy the prompt from `multiqc_data/llms-full.txt` into clipboard with the **Copy prompt** button in the toolbox, and use it with external services.

## Custom providers

Besides the built-in providers, MultiQC can use self-hosted models, other providers, or an internal gateway.
Set `ai_provider: custom` and either point MultiQC to an endpoint, or give it your own Python class.

### OpenAI-compatible endpoints

Most self-hosted servers (vLLM, Ollama, LM Studio, LiteLLM) and many providers implement the OpenAI chat completions API:

```yaml
ai_provider: custom
ai_model: your-model-name
ai_custom_endpoint: http://localhost:8000/v1/chat/completions
ai_custom_context_window: 32000 # Optional, see "Context window" above
```

Set the `OPENAI_API_KEY` environment variable to the endpoint's API key:

```bash
export OPENAI_API_KEY=your-api-key
```

To send the key in an `api-key` header instead of `Authorization: Bearer`, set `ai_auth_type: api-key`.

### Claude-compatible endpoints

For a gateway that speaks the Claude messages API, set `ai_custom_api_style: anthropic`,
and the API key in `ANTHROPIC_API_KEY`:

```yaml
ai_provider: custom
ai_custom_api_style: anthropic
ai_model: claude-opus-5-5
ai_custom_endpoint: https://llm-gateway.example.com/v1/messages
```

### Headers, request fields and responses

For endpoints with their own authentication or request and response formats:

```yaml
ai_custom_headers:
  x-api-key: ${MY_GATEWAY_TOKEN} # read from an environment variable
  Authorization: null # null removes a header MultiQC would otherwise send
ai_extra_query_options:
  model: null # null removes a request field, here the model name
  seed: 42 # any other field is added to the request
ai_custom_response_path: completion.content # where the answer is in the JSON response
```

When `ai_custom_headers` is set, the API key environment variable is optional.
Use `${VAR}` for secrets, so they stay out of config files.

### Your own client class

For anything else, write a Python class that subclasses `multiqc.core.ai.Client` and implements `query()`:

```python
# my_lab/llm.py
import os

import requests

from multiqc.core.ai import Client, LlmResponse


class GatewayClient(Client):
    name = "gateway"
    title = "Lab LLM gateway"
    model = "lab-model"

    def query(self, prompt: str) -> LlmResponse:
        response = requests.post(
            "https://llm.my-lab.org/generate",
            headers={"x-api-key": os.environ["LAB_LLM_TOKEN"]},
            json={"input": prompt},
            timeout=600,
        )
        response.raise_for_status()
        return LlmResponse(content=response.json()["output"], model=self.model)

    def context_window(self) -> int:
        return 200_000
```

```yaml
ai_provider: custom
ai_custom_client: my_lab.llm:GatewayClient
```

The class must be importable (installed, or on `PYTHONPATH`) and is created without arguments.
MultiQC still builds the prompt, fits it into `context_window()`, and anonymises sample names if configured.

### In the browser

In browser, you can also select "Custom" provider from the toolbox, and enter the endpoint URL, model name, API key,
optional extra query options, and context window manually.
In-browser summaries use the OpenAI chat completions format.
Custom headers and client classes are only used when generating summaries with MultiQC on the command line,
as they often contain secrets that should not be saved in a report.

## Configuring within Nextflow

If you're running MultiQC within a Nextflow pipeline, you probably don't want to edit the workflow code to configure AI summaries.
Most nf-core pipelines with MultiQC have a `--multiqc_config` option to provide an additional YAML config for MultiQC.
However, because API keys must be passed using environment variables, the recommended method is to use
environment vars for everything.

Using this approach means that no pipeline code needs adjustment, only a small addition to the Nextflow config.
For example, to use with OpenAI you would set the following in your Nextflow config:

```groovy
env {
  MULTIQC_AI_SUMMARY_FULL = 1     // Enable long summaries during report generation
  MULTIQC_AI_PROVIDER = "openai"  // Select OpenAI as provider
}
process.withName: MULTIQC {
  secret = [ 'OPEN_API_KEY' ]     // Access key for OpenAI
}
```

In this example, [Nextflow Secrets](https://nextflow.io/docs/latest/secrets.html) are used
to securely manage your API keys outside of your config file.
To add this Nextflow secret you would run the following command in the terminal:

```bash
$ nextflow secrets set OPENAI_API_KEY "xxxx"
```

Note that secrets can behave differently across different compute environment types.
See the [Nextflow docs](https://www.nextflow.io/docs/latest/secrets.html#process-directive)
for details.

:::tip

Save this Nextflow config to `~/.nextflow/config` and it
[will be applied](https://nextflow.io/docs/edge/config.html#configuration-file)
to every Nextflow pipeline you launch.

:::

## Security considerations

MultiQC AI summaries are used at your own risk.
Treat results with appropriate mistrust and consider what data you are sending to external services.

- API keys set in environment variables are not saved in report outputs
- API keys put in the toolbox are stored only in your browser's local storage
- No report data or keys are sent to any servers except the chosen AI provider

Seqera AI does not use inputs for subsequent fine-tuning or direct model improvement.
You can find our more information about Seqera's pledge for privacy at
[https://seqera.io/ai-trust/](https://seqera.io/ai-trust/)

### Sample anonymization

MultiQC provides an option to anonymize sample names when generating AI summaries, both during report generation and in the browser. This helps protect sensitive information when sharing summaries with AI providers.

To enable sample anonymization:

- For in-browser summaries: Toggle "Anonymize samples" in the AI toolbox section
- For report generation: Set `anonymize_samples: true` in your MultiQC config

When enabled, sample names are replaced with generic pseudonyms (e.g., "SAMPLE_1", "SAMPLE_2") before being sent to the AI provider. The anonymization is applied consistently across the entire report - each sample name gets the same pseudonym wherever it appears. When the AI response references samples, the pseudonyms are automatically converted back to the original sample names before displaying.

MultiQC replaces sample names that appear as typical keys in plots and tables, specifically:

- First column of a table table (e.g. the General statistics table)
- Labels of bars in bar plots
- Line names in line plots
- Data point labels in scatter plots
- Heatmap axis labels
- Violin plot data point names

But note that if a module creates some custom plot configuration where sample names are used elsewhere, anonymization would not be guaranteed.

:::info

Note that with the "Continue chat" button you would see the anonymized samples, which makes it less useful.

:::

[^seqera-ai-usage-limits]:
    Seqera Cloud Basic is free for small teams.
    It includes access to Seqera Co-Scientist, with a usage cap of 100 messages per calendar month.
    Researchers at qualifying academic institutions can apply for free access to Seqera Cloud Pro.
    See [Seqera Pricing](https://seqera.io/pricing/) for more details
