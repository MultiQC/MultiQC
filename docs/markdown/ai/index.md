---
title: AI summaries
description: Using AI to summarise MultiQC reports
---

# AI summaries

MultiQC v1.27 and newer can generate AI-powered summaries of your reports. These can be created at two points:

- When creating the report, baked into the report HTML.
- Dynamically in the browser, while viewing an existing HTML report.

The AI summaries are generated using LLMs (large-language models) AI, using
[Seqera AI](https://seqera.io/ask-ai/), [OpenAI](https://openai.com/) or [Anthropic](https://www.anthropic.com/).
MultiQC reports also have an option to copy a prompt to your clipboard, to paste into any provider you have access to.

:::warning

Never rely on AI-generated summaries. Whilst these summaries can be useful to get you started quickly with a report, they may give inaccurate analysis and miss important details.

:::

AI summaries work by sending report data to an LLM provider of your choice, via an API over the internet.
Be aware of what data you are sending, and to who.
For more information, see [Seqera AI: Your privacy matters](https://seqera.io/ai-trust/).

## Choosing a provider

To use native summary generation, MultiQC needs to communicate with an LLM provider's API.
All three supported services require an API key to work.
Remember: Treat your API keys like passwords and do not share them.

- [Seqera AI](https://seqera.io/ask-ai/)
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
- Any other provider, via your clipboard
  - You can use buttons in MultiQC reports to copy a prompt to your clipboard,
    in order to manually summarise report data.
    See [Copying prompts](#copying-prompts) for instructions.

Seqera AI is free to use.[^seqera-ai-usage-limits]
Use of OpenAI and Anthropic APIs are billed by their respective providers based on consumption.
Seqera AI uses the latest AI provider models under the hood, at the time of writing that is Anthropic Claude sonnet 3.5.

### Choosing a model

If you're using OpenAI or Anthropic you can choose the exact model used for report summaries.
This is done by setting `ai_model` in the MultiQC config.

- Anthropic model names must begin with `claude`
  - Default: `claude-3-5-sonnet-latest`.
  - Tested with Sonnet 3.5 and Haiku 3.5. See the [Anthropic docs](https://docs.anthropic.com/en/docs/intro-to-claude#model-options).
- OpenAI model names must being with `gpt`
  - Default: `gpt-4o`.
  - Tested with GPT-4o and GPT-4o-mini. See the [OpenAI docs](https://platform.openai.com/docs/models).

This model is used during report generation and also set as the default toolbox panel setting for browser report summaries.

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
  - `--ai-provider <provider>`: Choose AI provider. One of `seqera`, `openai` or `anthropic`. Default `seqera`
  - `--no-ai`: Disable AI toolbox and buttons in the report

- Alternatively, MultiQC configuration file:
  ```yaml
  ai_summary: false # Set to true for short summaries
  ai_summary_full: false # Set to true for  long summaries
  ai_provider: "seqera" # 'seqera', 'openai' or 'anthropic'. Default: 'seqera'
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
when processing input and generating responses, encompassing both the input prompt and the output.
At the time of writing, modern LLMs typically have a context window size in 128-200k tokens,
which translates to about 100-160k characters of report data.
That means that very large reports, of thousands of samples, might not fit in the available LLM context window.

MultiQC uses the following logic, moving on to the next step if the prompt is still too large:

1. Attempt to include all report data in the prompt.
2. Include just the general statistics table.
3. Include the general statistics table, without hidden-by-default columns.
4. Abort AI summary.

If you're unable to generate an AI summary, you can try the following:

- Hide additional columns in the general statistics table (see [Hiding Columns](../reports/customisation.md#hiding-columns)).
- Hide General statistics data in the browser, and request the AI summary dynamically:
  - Hide columns with the **Configure columns** button
  - Filter shown samples dynamically with the toolbox
- Copy the prompt from `multiqc_data/multiqc_ai_prompt.txt` into clipboard with the **Copy prompt** button in the toolbox, and use it with external services with a larger context window.

## Configuring within Nextflow

If you're running MultiQC within a Nextflow pipeline, you probably don't want to edit the workflow code to configure AI summaries.
Most nf-core pipelines with MultiQC have a `--multiqc_config` option to provide an additional YAML config for MultiQC.
However, because API keys must be passed using environment variables, the recommended method is to use
environment vars for everything.

Using this approach means that no pipeline code needs adjustment, only a small addition to the Nextflow config
by using the [`env` config scope](https://nextflow.io/docs/latest/reference/config.html#env).

For example, to use with OpenAI you would set the following in your Nextflow config:

```groovy
env {
   MULTIQC_AI_SUMMARY_FULL = 1              // Enable long summaries during report generation
   MULTIQC_AI_PROVIDER = "openai"           // Select OpenAI as provider
   OPENAI_API_KEY = secrets.OPENAI_API_KEY  // Access key for OpenAI
}
```

In this example, [Nextflow Secrets](https://nextflow.io/docs/latest/secrets.html) are used
to securely manage your API keys outside of your config file.
To add this Nextflow secret you would run the following command in the terminal:

```bash
$ nextflow secrets set OPENAI_API_KEY "xxxx"
```

:::tip

Save this Nextflow config to `~/.nextflow/config` and it
[will be applied](https://nextflow.io/docs/edge/config.html#configuration-file)
to every Nextflow pipeline you launch.

:::

### Using Seqera Platform

If using Seqera Platform, the Nextflow `env` config can be used when launching pipelines or adding them to the Launchpad.

Environment variables can also be added at the _Compute Environment_ level, which affects every pipeline
run using that CE without further modification.
This allows managing provider API keys at the workspace level.

To do this, toggle the **Environment variables** section when creating a Compute Environment and select **Add variable**.
Ensure that the **Target environment** has **Compute job** enabled.

![Seqera Platform: Setting environment variables at CE level](../../../docs/images/ai_seqera_platform_env_var_create.png)

Once added:

![Seqera Platform: Setting environment variables at CE level](../../../docs/images/ai_seqera_platform_env_vars.png)

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

MulitQC replaces sample names that appear as typical keys in plots and tables, specifically:

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
    It includes access to Seqera AI, with a usage cap of 100 messages per calendar month.
    Seqera AI usage is unlimited for Seqera Cloud Pro users.
    Researchers at qualifying academic institutions can apply for free access to Seqera Cloud Pro.
    See [Seqera Pricing](https://seqera.io/pricing/) for more details
