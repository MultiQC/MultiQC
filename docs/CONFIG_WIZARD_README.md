# MultiQC Configuration Wizard

The MultiQC Configuration Wizard is a self-contained HTML page that helps you build a `multiqc_config.yaml` file. Every option in `MultiQCConfig` is rendered as a form field with its description, default value, and inline YAML examples.

## How to use

Open the wizard file directly in any modern browser:

```bash
open docs/multiqc_config_wizard.html
```

No build step or local server is required. The page works fully offline once loaded; the only network dependencies are js-yaml, highlight.js, and Google Fonts.

Workflow:

1. Browse sections in the left sidebar, or type in the sidebar search.
2. Set the fields you want to override. Unset fields fall back to their defaults.
3. Click **Preview YAML** to see the output, or **Download Config** to save it as `multiqc_config.yaml`.
4. Drop the file in your project directory and run MultiQC as usual.

## Features

- **Per-field controls** chosen from the schema type:
  - Booleans and small enums (≤4 options) render as a segmented toggle.
  - Larger enums render as a dropdown.
  - `List[Literal[...]]` fields (eg. `export_plot_formats`) render as checkboxes.
  - Lists and objects render as a textarea pre-filled with the default value as YAML.
- **Inline examples**: fields that declare `examples=[...]` in their Pydantic `Field()` get a "Show example" link in the top-right of the form-group.
- **"Show advanced options" toggle** below the section nav reveals fields tagged in `UNCOMMON_PROPERTIES`. Searching always ignores this toggle, so any matching field is reachable.
- **Reset button (↻)** next to fields with pre-filled defaults; restores the original value.
- **Preview YAML mode** hides the editor and shows the generated YAML with syntax highlighting. Use **Editor** in the sidebar to go back.
- **Sticky beforeunload warning** when there are unsaved changes.
- **Empty-list semantics**: clearing a list field whose default is non-empty produces `field: []` in the YAML (an explicit override), not silent "unset".
- **Footer and header** record the MultiQC version used to build the wizard and the generation date.

## Architecture

- `multiqc/utils/config_schema.py` — Pydantic `MultiQCConfig` model. Source of truth for descriptions, types, defaults, examples, and `Literal` enums.
- `scripts/generate_config_wizard.py` — reads the schema and `multiqc/config_defaults.yaml`, decides which section each field belongs to, and substitutes the data into the template.
- `scripts/wizard_template.html` — the single-file HTML/CSS/JS template with `__MULTIQC_LOGO_SVG__`, `__CONFIG_DATA_JSON__`, `__MULTIQC_VERSION__`, and `__GENERATED_ON__` placeholders.
- `scripts/_config_schema_loader.py` — shared loader used by both the wizard and `generate_config_docs.py`.
- `docs/multiqc_config_wizard.html` — the rendered output (checked into the repo).

## Regenerating

```bash
python scripts/generate_config_wizard.py
```

The generator raises a `RuntimeError` if any schema property is missing from the wizard's `sections` map (or stale), so a new field can't slip through unnoticed. Run the tests after editing the schema:

```bash
pytest tests/test_config_wizard.py
```

Four tests guard the schema/wizard contract:

- every schema property is rendered (or explicitly in `SKIP_PROPERTIES`);
- `SKIP_PROPERTIES` and `UNCOMMON_PROPERTIES` only reference real fields;
- every key in `config_defaults.yaml` has a matching Pydantic field.

## Adding a new config option

1. Add the field to `MultiQCConfig` in `multiqc/utils/config_schema.py` with a `description=` and, where useful, `examples=[...]`. Use a `Literal[...]` (or `List[Literal[...]]`) when the value is restricted.
2. Add a default in `multiqc/config_defaults.yaml` if appropriate.
3. List the property in the matching section of `sections` in `scripts/generate_config_wizard.py`. Add it to `UNCOMMON_PROPERTIES` if it's internal or rarely set.
4. Regenerate the three artifacts:

   ```bash
   python scripts/generate_config_wizard.py
   python scripts/generate_config_docs.py
   python scripts/generate_config_schema.py
   ```

5. Run `pytest tests/test_config_wizard.py` to confirm the regression guards still pass.
