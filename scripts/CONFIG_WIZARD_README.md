# Config Wizard — contributor notes

End-user documentation lives at [`docs/markdown/getting_started/config_wizard.md`](../docs/markdown/getting_started/config_wizard.md). This file is for people working on the wizard's source.

## Files

- `multiqc/utils/config_schema.py` — Pydantic `MultiQCConfig` model. Source of truth for descriptions, types, defaults, examples, `Literal` enums, and per-field `section` / `uncommon` tags.
- `scripts/generate_config_wizard.py` — reads the schema and `multiqc/config_defaults.yaml`, classifies each field into a section, and substitutes the data into the template.
- `scripts/wizard_template.html` — the single-file HTML/CSS/JS template with `__MULTIQC_LOGO_SVG__`, `__CONFIG_DATA_JSON__`, `__CONFIG_SCHEMA_JSON__`, `__MULTIQC_VERSION__`, and `__GENERATED_ON__` placeholders.
- `scripts/_config_schema_loader.py` — shared loader used by both the wizard and `generate_config_docs.py`.
- `scripts/social_card.html` — standalone HTML used to render `docs/multiqc_config_wizard_social_card.png` via a headless screenshot at 1200×630.
- `docs/multiqc_config_wizard.html` — the rendered wizard (checked in).
- `docs/multiqc_config_wizard_social_card.png` — the Open Graph / Twitter Card image referenced from the rendered HTML's meta tags (also checked in).

## Regenerating

```bash
python scripts/generate_config_wizard.py
```

Run the tests after editing the schema:

```bash
pytest tests/test_config_wizard.py
```

Six tests guard the schema/wizard contract:

- every schema property declares a `section` (or is in `SKIP_PROPERTIES`);
- `SKIP_PROPERTIES` references only real schema fields;
- every user-facing field in `multiqc/config.py` is in the schema;
- `RUNTIME_CONFIG_ATTRS` (the allow-list of runtime-only attributes) is current;
- `config.py` type annotations agree with `MultiQCConfig` at the type-kind level;
- every key in `config_defaults.yaml` has a matching Pydantic field.

## Regenerating the social card

```bash
agent-browser open "file://$PWD/scripts/social_card.html"
agent-browser set viewport 1200 630
agent-browser screenshot docs/multiqc_config_wizard_social_card.png
```

Any tool that can take a 1200×630 viewport screenshot of a local HTML file works. The PNG and the wizard HTML are uploaded together to `seqera.io/multiqc_config_wizard` and `seqera.io/multiqc_config_wizard_social_card.png`.

## Adding a new config option

1. Add the field to `MultiQCConfig` in `multiqc/utils/config_schema.py` with `cfg(..., section="...", description=..., examples=[...])`. Use a `Literal[...]` (or `List[Literal[...]]`) when the value is restricted. Pass `uncommon=True` for fields that should hide behind the advanced toggle.
2. Add a default in `multiqc/config_defaults.yaml` if appropriate, and mirror the annotation in `multiqc/config.py`.
3. Regenerate the artifacts:

   ```bash
   python scripts/generate_config_wizard.py
   python scripts/generate_config_docs.py
   python scripts/generate_config_schema.py
   ```

4. Run `pytest tests/test_config_wizard.py`.
