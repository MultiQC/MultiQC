#!/usr/bin/env python3
"""
Generate markdown documentation from the MultiQC configuration Pydantic schema.

Reads the MultiQCConfig model and config_defaults.yaml, then produces a
Markdown reference grouped by logical section.  Run from the repo root::

    python scripts/generate_config_docs.py

Output: docs/markdown/config_schema.md
"""

import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union, get_args, get_origin, get_type_hints

import yaml


class _PrettierYamlDumper(yaml.SafeDumper):
    """YAML dumper that matches Prettier's output style.

    Two tweaks over the default ``SafeDumper``:

    1. Nested list dashes are indented one level deeper than the parent key,
       not at the parent's indent (``key:\\n  - item`` vs ``key:\\n- item``).
    2. Scalars that need quoting use double quotes instead of single. The
       emitter decides quoting style at write time, so we intercept it via
       ``choose_scalar_style`` rather than via a representer.
    """

    def increase_indent(self, flow=False, indentless=False):  # type: ignore[override]
        return super().increase_indent(flow, False)

    def choose_scalar_style(self):  # type: ignore[override]
        style = super().choose_scalar_style()
        return '"' if style == "'" else style


def _dump_yaml(value):
    """Dump YAML in a shape that round-trips through Prettier without changes."""
    return yaml.dump(
        value,
        Dumper=_PrettierYamlDumper,
        default_flow_style=False,
        sort_keys=False,
        allow_unicode=True,
        width=80,
    )


sys.path.insert(0, str(Path(__file__).parent))

from _config_schema_loader import load_schema_and_defaults, load_sections  # noqa: E402
from multiqc.utils.config_schema import (  # noqa: E402
    AiProviderLiteral,
    CleanPattern,
    GeneralStatsColumnConfig,
    GeneralStatsModuleConfig,
    MultiQCConfig,
    SearchPattern,
)

# Properties skipped in the section walk; sp is documented manually below under
# "Special Types".
SKIP_PROPERTIES = {"sp"}


def format_type_annotation(annotation):
    """Format a type annotation to a readable string."""
    if annotation is None:
        return "any"

    origin = get_origin(annotation)
    args = get_args(annotation)

    if origin is Union:
        # Drop the None branch; nearly every config field is Optional, so the
        # Optional[...] wrapper is noise in the rendered docs.
        non_null = [arg for arg in args if arg is not type(None)]
        if len(non_null) == 1:
            return format_type_annotation(non_null[0])
        formatted_args = [format_type_annotation(arg) for arg in non_null]
        return f"Union[{', '.join(formatted_args)}]"

    elif origin is list or origin is List:
        if args:
            return f"List[{format_type_annotation(args[0])}]"
        return "List"

    elif origin is dict or origin is Dict:
        if len(args) == 2:
            return f"Dict[{format_type_annotation(args[0])}, {format_type_annotation(args[1])}]"
        return "Dict"

    elif origin is Literal:
        literals = [f'"{a}"' if isinstance(a, str) else str(a) for a in args]
        return f"Literal[{', '.join(literals)}]"

    elif annotation is AiProviderLiteral:
        literals = [f"'{a}'" for a in get_args(annotation)]
        return f"Literal[{', '.join(literals)}]"

    # Handle basic types
    if annotation is str:
        return "str"
    elif annotation is int:
        return "int"
    elif annotation is float:
        return "float"
    elif annotation is bool:
        return "bool"
    elif annotation is dict or annotation is Dict:
        return "Dict"
    elif annotation is list or annotation is List:
        return "List"
    elif annotation is SearchPattern:
        return "SearchPattern"
    elif annotation is CleanPattern:
        return "CleanPattern"
    elif annotation is GeneralStatsColumnConfig:
        return "GeneralStatsColumnConfig"
    elif annotation is GeneralStatsModuleConfig:
        return "GeneralStatsModuleConfig"
    elif annotation is Any:
        return "Any"

    # For any other type, return its string representation
    return str(annotation).replace("typing.", "")


def format_default_value(value):
    """Format a default value for display in markdown."""
    if value is None:
        return "`None`"
    elif isinstance(value, bool):
        return f"`{str(value).lower()}`"
    elif isinstance(value, (int, float)):
        return f"`{value}`"
    elif isinstance(value, str):
        return f'`"{value}"`'
    elif isinstance(value, (list, dict)):
        # For complex objects, use JSON representation with some indentation
        # but keep it on one line for better markdown display
        return f"`{json.dumps(value, separators=(',', ':'))}`"
    else:
        return f"`{value}`"


# Inline-default character budget. Above this, the default is rendered
# below the type line as a collapsible <details> block so things like
# module_order's 200-entry list don't blow up the type line.
DEFAULT_INLINE_MAX = 120


def _is_uninformative_default(value):
    """True for defaults that add no information to the docs.

    Covers ``None``, empty lists/dicts, and dicts whose every value is itself
    uninformative (e.g. ``custom_content: {order: []}``).
    """
    if value is None:
        return True
    if isinstance(value, list) and not value:
        return True
    if isinstance(value, dict):
        if not value:
            return True
        return all(_is_uninformative_default(v) for v in value.values())
    return False


def render_default(value):
    """Return ``(inline, block)`` for use after a Type line.

    Uninformative defaults render as ``("", "")``. Short defaults go inline
    after Type. Long lists/dicts go into a collapsible block underneath so the
    Type line stays readable.
    """
    if _is_uninformative_default(value):
        return "", ""
    formatted = format_default_value(value)
    if len(formatted) <= DEFAULT_INLINE_MAX:
        return f" (default: {formatted})", ""
    yaml_text = _dump_yaml(value).rstrip()
    block = f"\n<details><summary>Default value</summary>\n\n```yaml\n{yaml_text}\n```\n\n</details>\n"
    return "", block


def generate_markdown_from_schema():
    """Generate markdown documentation from the MultiQC config schema."""
    config_attrs = dict(MultiQCConfig.__annotations__)
    properties, config_defaults, schema = load_schema_and_defaults()

    output = ["# MultiQC Configuration Reference\n\n"]
    output.append("""
This document describes all configuration options available in MultiQC.

## Introduction

MultiQC configuration can be set in several ways:

1. **Command line parameters** - Command line flags are available for many options (run `multiqc --help` to see all available options)
2. **Configuration files** - MultiQC looks for configuration files in the following locations (in order of precedence):
   - `<current working directory>/multiqc_config.yaml`
   - `~/.multiqc_config.yaml`
   - `<installation_dir>/multiqc/utils/config_defaults.yaml`
3. **Environment variables** - MultiQC checks for environment variables that match configuration options prefixed with `MULTIQC_`, for example: `MULTIQC_TITLE="My Report"`

Configuration values are loaded in the following order of precedence (highest to lowest):

1. Command line parameters
2. Current working directory config file
3. User home directory config file
4. Environment variables
5. Default configuration values

The options below can be specified in your YAML configuration files.
For boolean options, use `true` or `false` (all lowercase) in your YAML files.
""")

    # Group properties into logical sections from the schema's per-field tags.
    sections = load_sections(properties, skip=SKIP_PROPERTIES)

    # Generate markdown for each section
    for section, props in sections.items():
        if not props:
            continue

        output.append(f"## {section}\n")

        for prop_name in props:
            if prop_name in properties:
                prop = properties[prop_name]

                # Get type information from Pydantic model if available
                if prop_name in config_attrs:
                    type_info = format_type_annotation(config_attrs[prop_name])
                else:
                    # Fallback to JSON schema type
                    type_info = prop.get("type", "any")
                    if type_info == "array":
                        items = prop.get("items", {})
                        item_type = items.get("type", "any")
                        if "oneOf" in items:
                            item_types = [t.get("type", "any") for t in items.get("oneOf", [])]
                            item_type = " | ".join(item_types)
                        type_info = f"List[{item_type}]"
                    elif not type_info:
                        type_info = "any"

                # Get description
                description = prop.get("description", "")

                # Get default value - first from config_defaults.yaml, then fall back to schema
                if prop_name in config_defaults:
                    default_val = config_defaults[prop_name]
                else:
                    default_val = prop.get("default")
                default_inline, default_block = render_default(default_val)

                # Format the markdown
                output.append(f"### {prop_name}\n")
                output.append(f"**Type**: `{type_info}`{default_inline}\n")
                output.append(f"{description}\n")
                if default_block:
                    output.append(default_block)

                # Render examples as YAML code blocks
                examples = prop.get("examples") or []
                if examples:
                    label = "Example" if len(examples) == 1 else "Examples"
                    output.append(f"**{label}**:\n")
                    for ex in examples:
                        yaml_text = _dump_yaml({prop_name: ex}).rstrip()
                        output.append(f"```yaml\n{yaml_text}\n```\n")

        output.append("")  # Add blank line between sections

    # Describe special types
    output.append("## Special Types\n")

    # Search Pattern
    output.append("### SearchPattern\n")
    output.append("""Configuration for file search patterns used to find tool outputs.

The `SearchPattern` type is used in the `sp` configuration option to define patterns for finding and parsing tool output files.

Example:

```yaml
sp:
  fastqc:
    fn: "*_fastqc.zip"
  custom_tool:
    fn: "*.log"
    contents: "Started analysis"
```

Properties:\n\n""")
    search_pattern_attrs = get_type_hints(SearchPattern)
    for prop_name, prop_type in sorted(search_pattern_attrs.items()):
        # Get description from schema
        description = ""
        if "$defs" in schema and "SearchPattern" in schema["$defs"]:
            sp_props = schema["$defs"]["SearchPattern"].get("properties", {})
            if prop_name in sp_props:
                description = sp_props[prop_name].get("description", "")

        type_info = format_type_annotation(prop_type)
        output.append(f"- **{prop_name}** (`{type_info}`): {description}")
    output.append("")

    # Clean Pattern
    output.append("### CleanPattern\n")
    output.append("""Pattern for cleaning sample names.

The `CleanPattern` type is used in the `fn_clean_exts` and `extra_fn_clean_exts` configuration options to define patterns for cleaning sample names.

Example:

```yaml
fn_clean_exts:
  - type: truncate
    pattern: '_S\\d+_L\\d+'
  - type: regex
    pattern: '\\d{4}-\\d{2}-\\d{2}'
```

Properties:\n\n""")
    clean_pattern_attrs = get_type_hints(CleanPattern)
    for prop_name, prop_type in sorted(clean_pattern_attrs.items()):
        # Get description from schema
        description = ""
        if "$defs" in schema and "CleanPattern" in schema["$defs"]:
            cp_props = schema["$defs"]["CleanPattern"].get("properties", {})
            if prop_name in cp_props:
                description = cp_props[prop_name].get("description", "")

        type_info = format_type_annotation(prop_type)
        output.append(f"- **{prop_name}** (`{type_info}`): {description}")
    output.append("")

    # General Stats Column Configuration
    output.append("### GeneralStatsColumnConfig\n")
    output.append("""Configuration for columns in the general statistics table.

The `GeneralStatsColumnConfig` type is used in the `general_stats_columns` configuration option to customize the appearance and behavior of columns in the general statistics table.

Example:

```yaml
general_stats_columns:
  fastqc:
    columns:
      percent_duplicates:
        title: "% Dups"
        description: "Percentage of duplicate reads"
        scale: "RdYlGn-rev"
        max: 100
        min: 0
```

Properties:\n\n""")
    gs_col_attrs = get_type_hints(GeneralStatsColumnConfig)
    for prop_name, prop_type in sorted(gs_col_attrs.items()):
        # Get description from schema
        description = ""
        if "$defs" in schema and "GeneralStatsColumnConfig" in schema["$defs"]:
            gs_props = schema["$defs"]["GeneralStatsColumnConfig"].get("properties", {})
            if prop_name in gs_props:
                description = gs_props[prop_name].get("description", "")

        type_info = format_type_annotation(prop_type)
        output.append(f"- **{prop_name}** (`{type_info}`): {description}")

    text = "\n".join(output)
    # Match Prettier's expectations so the file survives the commit hook
    # without further edits: no trailing whitespace, at most one blank line
    # between blocks, single trailing newline.
    text = re.sub(r"[ \t]+(?=\n)", "", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.rstrip() + "\n"


if __name__ == "__main__":
    markdown = generate_markdown_from_schema()

    output_path = Path(__file__).parent.parent / "docs" / "markdown" / "config_schema.md"

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(markdown)
    except OSError as e:
        print(f"Error writing to {output_path}: {e}")
        sys.exit(1)

    print(f"Configuration documentation generated at {output_path}")
