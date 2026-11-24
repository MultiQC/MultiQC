---
title: Config Wizard
description: Build a multiqc_config.yaml in your browser, with live schema validation and a YAML editor.
---

# MultiQC Configuration Wizard

The MultiQC Configuration Wizard is a visual editor for `multiqc_config.yaml`. Every option in MultiQC appears as a form field with a description and inline examples. A live YAML editor next to the form stays in sync as you edit either side.

:::tip[Open the Config Wizard]

**[seqera.io/multiqc_config_wizard ↗](https://seqera.io/multiqc_config_wizard)**

:::

## What it does

- Live two-way sync
  - Form changes merge into the editor text, preserving comments and key order. Editor changes parse the YAML and push values back into the form once it's valid.
- Schema validation
  - Every option is checked against the same `MultiQCConfig` Pydantic schema MultiQC uses to load configs. Type mismatches and enum violations show as red squiggles in the editor; unknown keys get an amber squiggle with a "did you mean?" suggestion. Form rows pick up a matching coloured status label.
- Status filters
  - Chips at the top of the form help you to narrow it to just what's set, broken, unfilled, or matching the default.
- Hover docs
  - Hover over a key in the editor to see its description, type, and default.
- Click to reveal
  - Click a form row to highlight its line in the editor. Click a key in the editor to scroll the matching form row into view.

When you're done, hit **Copy** to put the YAML on the clipboard, or **Download** to save it as `multiqc_config.yaml`. Drop the file in your project directory (or anywhere MultiQC looks; see [Configuring MultiQC](config.md) for the search paths) and run MultiQC as usual.

## How to use it

1. Open [seqera.io/multiqc_config_wizard](https://seqera.io/multiqc_config_wizard)
2. Browse the sections in the left sidebar, search by key name, or paste an existing YAML into the editor to start from there.
3. Edit the fields you care about. Unset fields fall back to the MultiQC defaults; you don't need to touch them.
4. Copy or download the result.

You can leave the page and come back later: your in-progress YAML and form state persist in your browser's `localStorage`.

## Run it anywhere

The wizard is a single self-contained HTML file. It loads Monaco, Ajv, js-yaml, and Google Fonts from CDNs on first load; everything else lives in the page.

For offline use, on an air-gapped network, or pinned to your installed MultiQC version, download `docs/multiqc_config_wizard.html` from the [MultiQC repository](https://github.com/MultiQC/MultiQC) and open it in any modern browser:

```bash
open multiqc_config_wizard.html
```

The local file and the one on https://seqera.io come from the same template, regenerated on every MultiQC release.

## How it works

The wizard is generated from `MultiQCConfig`, the Pydantic model in `multiqc/utils/config_schema.py` that MultiQC itself uses at runtime. Each field's description, type, default, examples, and `Literal` enum flow through to the form, so the wizard tracks MultiQC release by release.

At build time, `scripts/generate_config_wizard.py` reads the schema, classifies each field into a section, and substitutes the resulting JSON into a single HTML template (`scripts/wizard_template.html`). The output is `docs/multiqc_config_wizard.html`, checked in alongside the source.

In the browser, [js-yaml](https://github.com/nodeca/js-yaml) parses the editor text into a JavaScript object on every change, and [Ajv](https://ajv.js.org/) validates that object against the JSON Schema exported from `MultiQCConfig`. Errors land in [Monaco's](https://microsoft.github.io/monaco-editor/) marker API as squiggles and hover tooltips.

The form ↔ editor sync is two-way. Form changes regenerate the YAML and merge them into the existing editor text, preserving comments where possible. Editor changes parse the YAML and push values back into the form widgets. Cycle guards stop the two from re-triggering each other.

The full source lives in the [MultiQC repository](https://github.com/MultiQC/MultiQC/tree/main/scripts).
