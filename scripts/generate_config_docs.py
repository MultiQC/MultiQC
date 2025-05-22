#!/usr/bin/env python3
"""
Script to generate markdown documentation from the MultiQC configuration schema.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union, get_args, get_origin, get_type_hints

import yaml

# Add parent directory to path so we can import multiqc
sys.path.insert(0, str(Path(__file__).parent.parent))

from multiqc.utils.config_schema import (
    AiProviderLiteral,
    CleanPattern,
    GeneralStatsColumnConfig,
    MultiQCConfig,
    SearchPattern,
)


def format_type_annotation(annotation):
    """Format a type annotation to a readable string."""
    if annotation is None:
        return "any"

    origin = get_origin(annotation)
    args = get_args(annotation)

    if origin is Union:
        # Handle Optional[Type] -> Union[Type, None]
        if type(None) in args:
            other_args = [arg for arg in args if arg is not type(None)]
            if len(other_args) == 1:
                return f"Optional[{format_type_annotation(other_args[0])}]"
            else:
                formatted_args = [format_type_annotation(arg) for arg in other_args]
                return f"Optional[Union[{', '.join(formatted_args)}]]"
        else:
            formatted_args = [format_type_annotation(arg) for arg in args]
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
        return "Literal['seqera', 'openai', 'anthropic', 'custom']"

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


def generate_markdown_from_schema():
    """Generate markdown documentation from the MultiQC config schema."""
    # Get Pydantic model attributes and type hints
    config_attrs = {name: field for name, field in MultiQCConfig.__annotations__.items()}

    # Get JSON schema for descriptions and default values
    schema = MultiQCConfig.model_json_schema()
    properties = schema.get("properties", {})

    # Load default values from config_defaults.yaml
    config_defaults_path = Path(__file__).parent.parent / "multiqc" / "config_defaults.yaml"
    with open(config_defaults_path, "r") as f:
        config_defaults = yaml.safe_load(f)

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

    # Group properties into logical sections
    sections = {
        "General": [],
        "Report Appearance": [
            "title",
            "subtitle",
            "intro_text",
            "report_comment",
            "report_header_info",
            "show_analysis_paths",
            "show_analysis_time",
            "custom_logo",
            "custom_logo_url",
            "custom_logo_title",
            "custom_css_files",
            "simple_output",
            "template",
        ],
        "Output Options": [
            "output_fn_name",
            "data_dir_name",
            "plots_dir_name",
            "data_format",
            "force",
            "make_data_dir",
            "zip_data_dir",
            "data_dump_file",
            "data_dump_file_write_raw",
            "export_plots",
            "export_plots_timeout",
            "make_report",
            "make_pdf",
        ],
        "MegaQC Integration": ["megaqc_url", "megaqc_access_token", "megaqc_timeout"],
        "AI Summary": [
            "ai_summary",
            "ai_summary_full",
            "ai_provider",
            "ai_model",
            "ai_custom_endpoint",
            "ai_auth_type",
            "ai_retries",
            "ai_extra_query_options",
            "ai_custom_context_window",
            "ai_prompt_short",
            "ai_prompt_full",
            "no_ai",
            "ai_anonymize_samples",
        ],
        "Seqera Integration": ["seqera_api_url", "seqera_website"],
        "Plot Settings": [
            "plots_force_flat",
            "plots_force_interactive",
            "plots_export_font_scale",
            "plots_flat_numseries",
            "plots_defer_loading_numseries",
            "lineplot_number_of_points_to_hide_markers",
            "barplot_legend_on_bottom",
            "violin_downsample_after",
            "violin_min_threshold_outliers",
            "violin_min_threshold_no_points",
        ],
        "Table Settings": [
            "collapse_tables",
            "max_table_rows",
            "max_configurable_table_columns",
            "decimalPoint_format",
            "thousandsSep_format",
        ],
        "Sample Names": [
            "prepend_dirs",
            "prepend_dirs_depth",
            "prepend_dirs_sep",
            "fn_clean_sample_names",
            "use_filename_as_sample_name",
            "fn_clean_exts",
            "fn_clean_trim",
            "extra_fn_clean_exts",
            "extra_fn_clean_trim",
            "sample_names_ignore",
            "sample_names_ignore_re",
            "sample_names_only_include",
            "sample_names_only_include_re",
            "sample_names_rename_buttons",
            "sample_names_replace",
            "sample_names_replace_regex",
            "sample_names_replace_exact",
            "sample_names_replace_complete",
            "sample_names_rename",
        ],
        "Toolbox": [
            "show_hide_buttons",
            "show_hide_patterns",
            "show_hide_regex",
            "show_hide_mode",
            "highlight_patterns",
            "highlight_colors",
            "highlight_regex",
        ],
        "Performance & Debugging": [
            "profile_runtime",
            "profile_memory",
            "verbose",
            "no_ansi",
            "quiet",
            "lint",
            "strict",
            "development",
            "log_filesize_limit",
            "filesearch_lines_limit",
            "report_readerrors",
        ],
        "File Discovery": [
            "require_logs",
            "ignore_symlinks",
            "ignore_images",
            "fn_ignore_dirs",
            "fn_ignore_paths",
            "filesearch_file_shared",
        ],
        "Other": [],
    }

    # Assign properties to sections
    for prop_name in properties:
        assigned = False
        for section, props in sections.items():
            if prop_name in props:
                assigned = True
                break
        if not assigned and prop_name not in ["$schema", "$defs"]:
            sections["Other"].append(prop_name)

    # Generate markdown for each section
    for section, props in sections.items():
        if not props:
            continue

        output.append(f"## {section}\n")

        for prop_name in sorted(props):
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
                default = ""
                if prop_name in config_defaults:
                    default_val = config_defaults[prop_name]
                    default = f" (default: {format_default_value(default_val)})"
                elif "default" in prop:
                    default_val = prop["default"]
                    if default_val is None:
                        default = " (default: `None`)"
                    elif isinstance(default_val, (list, dict)):
                        default = f" (default: `{json.dumps(default_val)}`)"
                    else:
                        default = f" (default: `{default_val}`)"

                # Format the markdown
                output.append(f"### {prop_name}\n")
                output.append(f"**Type**: `{type_info}`{default}\n")
                output.append(f"{description}\n")

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
    fn: '*_fastqc.zip'
  custom_tool:
    fn: '*.log'
    contents: 'Started analysis'
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

    return "\n".join(output)


if __name__ == "__main__":
    # Generate markdown
    markdown = generate_markdown_from_schema()

    # Output file
    output_path = Path(__file__).parent.parent / "docs" / "markdown" / "config_schema.md"

    # Write to file
    with open(output_path, "w") as f:
        f.write(markdown)

    print(f"Configuration documentation generated at {output_path}")
