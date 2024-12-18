---
title: AI Summaries
description: Using AI to summarise MultiQC reports
---

# AI Summaries

MultiQC can generate AI-powered summaries of your reports. These can be created at two points:

- When creating the report (baked into the report HTML)
- Dynamically in the browser, whilst viewing an existing HTML report (ephemeral)

The AI summaries are generated using LLMs (large-language models) AI, using
[Seqera AI](https://seqera.io/ask-ai/), [OpenAI](https://openai.com/) or [Anthropic](https://www.anthropic.com/).
MultiQC reports also have an option to copy a prompt to your clipboard, to paste into any provider you have access to.

:::warning

Never rely on AI-generated summaries. Whilst these summaries can be useful to get you started quickly with a report, they may give inaccurate analysis and miss important details.

:::

:::warning

AI summaries work by sending report data to an LLM provider of your choice, via an API over the internet.
Be aware of what data you are sending, and to who.

:::

## Choosing a provider

To use native summary generation, MultiQC needs to communicate with an LLM provider's API.
All three supported services require an API key to work.

Seqera AI is free to use, though there are usage limits.
OpenAI / Anthropic API usage is billed based on consumption.

1. [Seqera AI](https://seqera.io/ask-ai/)
   - Register for free at [seqera.io](https://seqera.io/)
   - Create a new key on the _"User tokens"_ page: [https://cloud.seqera.io/tokens](https://cloud.seqera.io/tokens)
2. [OpenAI](https://openai.com/)
   - Register at [platform.openai.com](https://platform.openai.com/signup) (NB: different to ChatGPT)
   - Add a payment method to your account to enable API usage beyond any trial credits
   - Create a new secret key on the _API Keys_ section under your profile
3. [Anthropic](https://www.anthropic.com/)
   - Sign up at [https://console.anthropic.com](https://console.anthropic.com)
   - Add a payment method to enable API access
   - Create a new key on on the _API Keys_ section in your account settings

:::note

Treat your API key like a password and do not share it.

:::

### Choosing a model

If you're using OpenAI or Anthropic you can choose the exact model used for report summaries.
This is done by setting `ai_model` in the MultiQC config.

- Anthropic:
  - Model names must begin with `claude`
  - Default: `claude-3-5-sonnet-latest`. Tested with `claude-3-5-sonnet` and `claude-3-5-haiku`.
- OpenAI:
  - Model names must being with `gpt`
  - Default: `gpt-4o`. Tested with `gpt-4o`.

This model will then be used during report generation and also set as the default toolbox panel setting for browser report summaries.

## Summaries during report generation

MultiQC can generate AI summaries at run time, when generating reports.
Summary text will be included within the report HTML as static text and will be visible to anyone viewing the report,
however it is shared.

AI summaries are disabled by default when running MultiQC.
To use, you must enable them either on the command line or via a MultiQC config file.

1. Command line flags:

   - `--ai` / `--ai-summary`: Generate a short report summary and put it on top of the report (fast)
   - `--ai-summary-full`: Generate a detailed version of the summary with analysis and recommendations (slower)
   - `--ai-provider <provider>`: Choose AI provider. `seqera`, `openai` or `anthropic`. Default `seqera`.

2. Alternatively, MultiQC configuration file:
   ```yaml
   ai_summary: false # Set to 'true' for short summaries
   ai_summary_full: false # Set to 'true' for  long summaries
   ai_provider: "seqera" # 'seqera', 'openai' or 'anthropic'. Default: 'seqera'
   ```

You will need to find and set your provider's API key in order to access its service
_(see [Choosing a provider](#choosing-a-provider) for how to get an API key)_.
The API keys are supplied by setting the following environment variables:

```bash
export SEQERA_ACCESS_TOKEN="..."  # or TOWER_ACCESS_TOKEN
export OPENAI_API_KEY="..."
export ANTHROPIC_API_KEY="..."
```

Only one key is needed, according to which service you want to use.
If you run MultiQC without the appropriate key you will get a warning printed to the console,
but report generation will otherwise proceed without the summary and MultiQC will not return an error exit code.

MultiQC uses the [python-dotenv](https://saurabh-kumar.com/python-dotenv/) package,
so you can also use an `.env` file either in the current working directory or the MultiQC source code directory.

:::note

These environment variables will only be used for `--ai-summary`/`--ai-summary-full` generation.
They are not saved by MultiQC and cannot be used for summary generation within reports.

:::

## In-browser AI Summaries

In addition to summaries during report generation, MultiQC can also create summaries dynamically in reports.
This can be useful as the person viewing a report is often different to the person who generated it.
It also means that summaries can be generated on demand, only when needed.

Summaries generated when browsing reports are _ephemeral_ - that is, they are not saved in the HTML.
If you generate a summary and share the report then others will not see it.
MultiQC tries to save the summary response within your browser's [local storage](https://www.w3schools.com/html/html5_webstorage.asp)
so that it shows the next time you open the same report, but this process is imperfect and may not always work.

### Configuring the AI provider

When you first try to generate a summary in the browser, you must supply the LLM provider's API key.

Click either the the "Summarize report" button in the top right corner of the report, or the AI toolbox icon:

![ai_summarize_report](../../../docs/images/ai_empty.png)

Choose an AI provider and enter the relevant API key
_(see [Choosing a provider](#choosing-a-provider) for how to get an API key)_.

API keys will be stored only in your browser's local storage, so will not be shared with anyone if you send the HTML report to someone else.
They are used to send report data directly to your AI provider of choice.

![ai_summarize_report_toolbox_keys](../../../docs/images/ai_toolbox_keys.png)

### Summarising the report

Once your provider API key is configured, click _"Summarize report"_ to generate an overview summary of the entire report.

:::info

The summary text is interactive! Click an underlined sample name to highlight that sample throughout the report.

:::

![ai_summarize_report_toolbox](../../../docs/images/ai_toolbox.png)

### Section-level summaries

Besides a global report-level AI summary, you can generate a summary for each plot or table separately using buttons next to each section:

![ai_summarize_buttons](../../../docs/images/ai_summarize_buttons.png)

You can also copy the entire prompt that MultiQC would use to generate a summary. Using the toolbox, you can enable the "Copy prompt" buttons to be shown next to each "Summarize" button:

![ai_toolbox_copy_checkbox](../../../docs/images/ai_toolbox_copy_checkbox.png)

![ai_copy_buttons](../../../docs/images/ai_copy_buttons.png)

A button will copy the LLM-friendly formatted report data along with the system prompt into your clipboard, which you can then paste into your an AI chat interace of your provider of choice.

You can also copy the prompt for the entire report by clicking the "Copy prompt" button in the toolbox:

![ai_toolbox_copy_button](../../../docs/images/ai_toolbox_copy_button.png)

## Context window

At the time of writing, modern LLMs typically have a context window size in `128,000` - `200,000` tokens, which translates to about `100,000` - `160,000` characters in report data. That means that very large reports - of thousands of samples - might not fit the context window. If the entire report doesn't fit the context window, only the general statistics table will be included in the prompt, as it represents the essential overall information of the report. If it doesn't fit, hidden-by-default columns be excluded from the prompt. If even that doesn't help, it wouldn't be possible to generate an AI summary - however, you can try the following:

- Hide some columns in the general statistics table (see https://docs.seqera.io/multiqc/reports/customisation#hiding-columns) to reduce the context.
- Open the HTML report in the browser, hide columns with the "Configure columns" button, or filter shown samples dynamically with the toolbox, and request the AI summary dynamically.
- Copy the prompt from `multiqc_data/multiqc_ai_prompt.txt` or into clipboard with the "Copy prompt" button in the toolbox, and use it with extrenal services with a larger context window.

## Security Considerations

- API keys put in the toolbox are stored only in your browser's local storage
- No report data or keys are sent to any servers except the chosen AI provider
- However, consider security implications when sending reports with sensitive data to your AI provider
