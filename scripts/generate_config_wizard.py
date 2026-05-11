#!/usr/bin/env python3
"""
Generate an interactive HTML configuration wizard for MultiQC.

Reads the MultiQCConfig Pydantic model and config_defaults.yaml, then
produces a single self-contained HTML page with a guided form for building
a ``multiqc_config.yaml`` file.  Run from the repo root::

    python scripts/generate_config_wizard.py

Output: docs/multiqc_config_wizard.html
"""

import json
import sys
from pathlib import Path
from typing import get_args

import yaml

# Add parent directory to path so we can import multiqc
sys.path.insert(0, str(Path(__file__).parent.parent))

from multiqc.utils.config_schema import AiProviderLiteral, MultiQCConfig

MULTIQC_LOGO_SVG = """\
<svg width="1318" height="250" viewBox="0 0 1318 250" fill="none" xmlns="http://www.w3.org/2000/svg">
<path d="M46.08 119.61C48.7 80.16 80.39 48.55 119.88 46.07V0C54.92 2.56 2.7 54.68 0 119.61H46.08Z" fill="#F18046"/>
<path d="M119.61 203.919C80.16 201.299 48.55 169.609 46.07 130.119H0C2.56 195.079 54.68 247.299 119.61 249.999V203.919Z" fill="#F18046"/>
<path d="M130.389 46.08C169.839 48.7 201.449 80.39 203.929 119.88H249.999C247.439 54.92 195.319 2.7 130.389 0V46.08Z" fill="#F18046"/>
<path d="M249.999 203.919C210.549 201.299 178.939 169.609 176.459 130.119H130.389C132.949 195.079 185.069 247.299 249.999 249.999V203.919Z" fill="#F18046"/>
<path d="M1100.98 249.999V179.569H1096.12C1093.89 183.629 1090.65 187.529 1086.39 191.269C1082.13 195.019 1076.71 198.109 1070.13 200.539C1063.54 202.969 1055.39 204.189 1045.66 204.189C1033.09 204.189 1021.54 201.149 1011 195.069C1000.46 188.989 992.05 180.229 985.77 168.769C979.48 157.319 976.35 143.489 976.35 127.269V122.709C976.35 106.499 979.54 92.6691 985.93 81.2091C992.31 69.7591 1000.77 60.9891 1011.31 54.9091C1021.85 48.8291 1033.3 45.7891 1045.66 45.7891C1060.25 45.7891 1071.45 48.4291 1079.25 53.6891C1087.05 58.9591 1092.88 64.9391 1096.73 71.6291H1101.59V45.7891H1132.29V249.979H1100.98V249.999ZM1054.47 176.829C1068.25 176.829 1079.5 172.469 1088.21 163.759C1096.92 155.049 1101.28 142.579 1101.28 126.369V123.629C1101.28 107.619 1096.87 95.2591 1088.06 86.5391C1079.24 77.8291 1068.04 73.4691 1054.47 73.4691C1040.9 73.4691 1030 77.8291 1021.18 86.5391C1012.36 95.2591 1007.96 107.619 1007.96 123.629V126.369C1007.96 142.589 1012.37 155.049 1021.18 163.759C1030 172.479 1041.09 176.829 1054.47 176.829Z" fill="#333" class="mqc-logo-text"/>
<path d="M1235.96 204.189C1221.57 204.189 1208.55 201.149 1196.9 195.069C1185.24 188.989 1176.02 180.169 1169.24 168.619C1162.45 157.069 1159.06 143.189 1159.06 126.969V123.019C1159.06 106.809 1162.45 92.9793 1169.24 81.5193C1176.03 70.0693 1185.25 61.2593 1196.9 55.0693C1208.55 48.8893 1221.57 45.7993 1235.96 45.7993C1250.35 45.7993 1262.61 48.4393 1272.74 53.6993C1282.87 58.9693 1291.03 65.9593 1297.21 74.6793C1303.39 83.3993 1307.39 93.0193 1309.22 103.559L1278.82 109.939C1277.8 103.249 1275.68 97.1693 1272.44 91.6993C1269.2 86.2293 1264.64 81.8693 1258.76 78.6293C1252.88 75.3893 1245.48 73.7693 1236.57 73.7693C1227.66 73.7693 1220 75.7493 1213.01 79.6993C1206.02 83.6493 1200.49 89.3293 1196.44 96.7193C1192.38 104.119 1190.36 113.089 1190.36 123.619V126.359C1190.36 136.899 1192.38 145.919 1196.44 153.419C1200.49 160.919 1206.02 166.599 1213.01 170.439C1220 174.289 1227.85 176.219 1236.57 176.219C1249.74 176.219 1259.77 172.829 1266.67 166.039C1273.56 159.249 1277.92 150.589 1279.74 140.049L1310.14 147.039C1307.71 157.379 1303.4 166.909 1297.22 175.619C1291.04 184.339 1282.88 191.279 1272.75 196.439C1262.61 201.609 1250.35 204.189 1235.97 204.189H1235.96Z" fill="#333" class="mqc-logo-text"/>
<path d="M330.02 204.19V45.8096H360.72V65.8696H365.58C368.42 60.5996 372.98 55.9396 379.26 51.8896C385.54 47.8396 394.05 45.8096 404.8 45.8096C415.55 45.8096 424.91 48.0896 431.7 52.6496C438.49 57.2096 443.6 63.0396 447.05 70.1296H451.91C455.35 63.2396 460.37 57.4696 466.96 52.7996C473.54 48.1396 482.92 45.8096 495.08 45.8096C504.81 45.8096 513.42 47.7896 520.92 51.7396C528.42 55.6896 534.4 61.5696 538.86 69.3696C543.32 77.1696 545.55 86.8496 545.55 98.3996V204.19H514.24V100.83C514.24 91.7096 511.76 84.6696 506.79 79.6996C501.82 74.7396 494.78 72.2496 485.66 72.2496C475.93 72.2496 468.13 75.3896 462.25 81.6696C456.37 87.9496 453.43 96.9696 453.43 108.73V204.19H422.12V100.83C422.12 91.7096 419.64 84.6696 414.67 79.6996C409.7 74.7396 402.66 72.2496 393.54 72.2496C383.81 72.2496 376.01 75.3896 370.13 81.6696C364.25 87.9496 361.31 96.9696 361.31 108.73V204.19H330H330.02Z" fill="#333" class="mqc-logo-text"/>
<path d="M634.529 203.85C623.179 203.85 613.049 201.31 604.129 196.25C595.209 191.19 588.219 183.99 583.149 174.67C578.079 165.35 575.549 154.2 575.549 141.23V45.8096H606.859V139.1C606.859 152.07 610.099 161.65 616.589 167.83C623.069 174.01 632.089 177.1 643.649 177.1C656.419 177.1 666.699 172.8 674.509 164.18C682.309 155.57 686.209 143.16 686.209 126.94V45.8096H717.519V204.19H686.819V178.62H681.959C679.119 184.7 674.059 190.43 666.759 195.8C659.459 201.17 648.719 203.86 634.539 203.86L634.529 203.85Z" fill="#333" class="mqc-logo-text"/>
<path d="M747.521 204.19V0H778.831V204.19H747.521Z" fill="#333" class="mqc-logo-text"/>
<path d="M915.029 204.19V45.8096H946.339V204.19H915.029Z" fill="#333" class="mqc-logo-text"/>
<path d="M885.03 71.47V45.81H865.56L848.45 62.92H840.26V0H808.83V173.98C808.83 183.13 811.52 190.46 816.92 195.95C822.31 201.44 829.67 204.19 839.04 204.19H885.02V177.95H849.1C843.2 177.95 840.25 174.9 840.25 168.8V71.48H885.01L885.03 71.47Z" fill="#333" class="mqc-logo-text"/>
</svg>"""


def generate_config_wizard():
    """Generate a self-contained HTML configuration wizard for MultiQC.

    The wizard shows all options on a single scrollable page, with a sidebar
    for navigation and a global search bar.  Boolean fields use a tri-state
    select (not set / true / false) so that ``default: true`` values are
    handled clearly.
    """

    schema = MultiQCConfig.model_json_schema()
    properties = schema.get("properties", {})

    config_defaults_path = Path(__file__).parent.parent / "multiqc" / "config_defaults.yaml"
    try:
        with open(config_defaults_path, "r") as f:
            config_defaults = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Could not find {config_defaults_path}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML: {e}")
        sys.exit(1)

    sections = {
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
        "Sample Names": [
            "prepend_dirs",
            "prepend_dirs_depth",
            "prepend_dirs_sep",
            "fn_clean_sample_names",
            "use_filename_as_sample_name",
            "sample_names_ignore",
            "sample_names_ignore_re",
            "sample_names_only_include",
            "sample_names_only_include_re",
        ],
        "File Discovery": [
            "require_logs",
            "ignore_symlinks",
            "ignore_images",
            "fn_ignore_dirs",
            "fn_ignore_paths",
        ],
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
    }

    config_data: dict = {}
    for section_name, section_props in sections.items():
        config_data[section_name] = {}
        for prop_name in section_props:
            if prop_name not in properties:
                continue
            prop = properties[prop_name]

            prop_type = prop.get("type", "string")
            if "anyOf" in prop:
                type_options = [t.get("type") for t in prop.get("anyOf", []) if "type" in t]
                if type_options:
                    prop_type = type_options[0]

            enum_values = None
            if "enum" in prop:
                enum_values = prop["enum"]
            elif prop_name == "ai_provider":
                enum_values = list(get_args(AiProviderLiteral))

            default_val = config_defaults.get(prop_name)

            config_data[section_name][prop_name] = {
                "type": prop_type,
                "description": prop.get("description", ""),
                "default": default_val,
                "enum": enum_values,
            }

    # Escape JSON data for safe embedding inside a <script> tag
    config_json = json.dumps(config_data, indent=8)
    config_json_escaped = config_json.replace("</", "<\\/")

    html_content = _build_html(config_json_escaped)
    return html_content


def _build_html(config_json_escaped: str) -> str:
    """Return the complete HTML string for the wizard."""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MultiQC Configuration Wizard</title>
    <link rel="icon" type="image/svg+xml" href="data:image/svg+xml;utf8,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 250 250'%3E%3Cpath fill='%23F18046' d='M46.08 119.61C48.7 80.16 80.39 48.55 119.88 46.07V0C54.92 2.56 2.7 54.68 0 119.61H46.08Z'/%3E%3Cpath fill='%23F18046' d='M119.61 203.919C80.16 201.299 48.55 169.609 46.07 130.119H0C2.56 195.079 54.68 247.299 119.61 249.999V203.919Z'/%3E%3Cpath fill='%23F18046' d='M130.389 46.08C169.839 48.7 201.449 80.39 203.929 119.88H249.999C247.439 54.92 195.319 2.7 130.389 0V46.08Z'/%3E%3Cpath fill='%23F18046' d='M249.999 203.919C210.549 201.299 178.939 169.609 176.459 130.119H130.389C132.949 195.079 185.069 247.299 249.999 249.999V203.919Z'/%3E%3C/svg%3E">
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,500;12..96,600;12..96,700&family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <style>
        :root {{
            /* Seqera tokens */
            --navy: #201637;            /* --color-brand */
            --navy-2: #2d273c;          /* --color-brand-1000 */
            --ink: #160f26;             /* heading ink */
            --ink-mute: #5c5767;        /* body text */
            --ink-soft: #736f7d;        /* caption */
            --ink-faint: #a29fa8;       /* subtle */
            --teal: #31c9ac;            /* --color-nextflow-500 */
            --teal-border: #0cae8e;     /* --color-nextflow-700 */
            --teal-border-active: #087f68; /* --color-nextflow-900 */
            --teal-tint: #e2f7f3;
            --orange: #f18046;          /* --color-multiqc-600 (logo only) */
            --paper: #ffffff;
            --surface-soft: #f8f8f8;    /* --color-brand-50 */
            --border: rgba(0, 0, 0, 0.15);
            --border-soft: #e8e7e9;     /* --color-brand-200 */
            --focus: #31c9ac;
        }}

        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}

        html, body {{
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            color: var(--ink-mute);
            background: var(--paper);
            font-size: 16px;
            line-height: 1.55;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
        }}

        .display {{
            font-family: 'Bricolage Grotesque', 'Inter', sans-serif;
            font-variation-settings: 'opsz' 32;
            color: var(--ink);
        }}

        code, pre, .mono {{
            font-family: 'JetBrains Mono', ui-monospace, SFMono-Regular, Monaco, monospace;
        }}

        a {{
            color: var(--navy);
            text-decoration: underline;
            text-underline-offset: 3px;
            text-decoration-thickness: 1px;
        }}
        a:hover {{ color: var(--teal-border); }}

        /* Fixed left side-nav */
        .sidebar {{
            position: fixed;
            top: 0;
            left: 0;
            width: 280px;
            height: 100vh;
            background: var(--navy);
            color: rgba(255, 255, 255, 0.85);
            display: flex;
            flex-direction: column;
            z-index: 50;
            border-right: 1px solid rgba(255, 255, 255, 0.06);
        }}
        .sidebar-brand {{
            padding: 22px 22px 18px;
            border-bottom: 1px solid rgba(255, 255, 255, 0.08);
        }}
        .sidebar-brand svg {{ height: 26px; width: auto; display: block; }}
        .sidebar-brand svg path.mqc-logo-text {{ fill: white; }}
        .sidebar-eyebrow {{
            margin-top: 12px;
            font-size: 0.7em;
            font-weight: 600;
            letter-spacing: 0.16em;
            text-transform: uppercase;
            color: var(--teal);
        }}

        .sidebar-search {{
            padding: 18px 18px 12px;
            position: relative;
        }}
        .sidebar-search input {{
            width: 100%;
            padding: 10px 14px 10px 38px;
            border: 1px solid rgba(255, 255, 255, 0.15);
            border-radius: 4px;
            background-color: rgba(255, 255, 255, 0.05);
            color: white;
            font-family: inherit;
            font-size: 0.9em;
            background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='14' height='14' viewBox='0 0 24 24' fill='none' stroke='%23A29FA8' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Ccircle cx='11' cy='11' r='7'/%3E%3Cpath d='M20.5 20.5l-4.4-4.4'/%3E%3C/svg%3E");
            background-repeat: no-repeat;
            background-position: 14px center;
            transition: border-color 0.15s ease, background-color 0.15s ease;
        }}
        .sidebar-search input::placeholder {{ color: rgba(255, 255, 255, 0.45); }}
        .sidebar-search input:focus {{
            outline: none;
            border-color: var(--teal);
            background-color: rgba(255, 255, 255, 0.08);
        }}

        .sidebar-sections {{
            flex: 1;
            overflow-y: auto;
            padding: 6px 12px 16px;
        }}
        .sidebar-sections h3 {{
            color: rgba(255, 255, 255, 0.45);
            margin: 8px 10px 8px;
            font-size: 0.7em;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.14em;
        }}
        .section-nav {{ list-style: none; }}
        .section-nav li {{ margin-bottom: 2px; }}
        .section-nav a {{
            color: rgba(255, 255, 255, 0.78);
            text-decoration: none;
            padding: 8px 12px;
            border-radius: 4px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-size: 0.92em;
            font-weight: 500;
            transition: background 0.15s ease, color 0.15s ease;
        }}
        .section-nav a:hover {{
            background: rgba(255, 255, 255, 0.06);
            color: white;
        }}
        .section-nav a.active {{
            background: rgba(49, 201, 172, 0.15);
            color: white;
            box-shadow: inset 2px 0 0 var(--teal);
        }}
        .section-nav .badge {{
            background: rgba(255, 255, 255, 0.1);
            color: rgba(255, 255, 255, 0.7);
            border-radius: 999px;
            padding: 2px 8px;
            font-size: 0.74em;
            font-weight: 500;
            min-width: 22px;
            text-align: center;
        }}
        .section-nav a.active .badge {{
            background: var(--teal);
            color: var(--navy);
        }}

        .sidebar-actions {{
            border-top: 1px solid rgba(255, 255, 255, 0.08);
            padding: 16px 18px 14px;
            display: flex;
            flex-direction: column;
            gap: 8px;
        }}
        .sidebar-actions .btn {{
            width: 100%;
            text-align: center;
        }}
        .sidebar-actions .btn-primary {{
            background: transparent;
            color: white;
            border-color: rgba(255, 255, 255, 0.3);
        }}
        .sidebar-actions .btn-primary:hover {{
            background: white;
            color: var(--navy);
            border-color: white;
        }}
        .sidebar-actions .btn-secondary {{
            background: transparent;
            color: rgba(255, 255, 255, 0.6);
            border-color: rgba(255, 255, 255, 0.12);
            padding: 6px 14px;
            font-size: 0.85em;
        }}
        .sidebar-actions .btn-secondary:hover {{
            background: rgba(255, 255, 255, 0.05);
            color: white;
            border-color: rgba(255, 255, 255, 0.25);
        }}

        .sidebar-links {{
            border-top: 1px solid rgba(255, 255, 255, 0.08);
            padding: 12px 18px 18px;
            display: flex;
            flex-direction: column;
            gap: 2px;
        }}
        .sidebar-links a {{
            color: rgba(255, 255, 255, 0.55);
            text-decoration: none;
            font-size: 0.82em;
            padding: 5px 10px;
            border-radius: 4px;
            transition: background 0.15s ease, color 0.15s ease;
        }}
        .sidebar-links a:hover {{
            background: rgba(255, 255, 255, 0.06);
            color: white;
        }}
        .sidebar-links a::before {{
            content: "↗";
            display: inline-block;
            margin-right: 8px;
            color: var(--teal);
            font-size: 0.85em;
        }}

        /* Right-hand main column */
        .main {{
            margin-left: 280px;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }}

        /* Hero in the main column */
        .header {{
            background: var(--paper);
            color: var(--ink);
            padding: 48px 56px 36px;
            border-bottom: 1px solid var(--border-soft);
        }}
        .header h1 {{
            font-family: 'Bricolage Grotesque', 'Inter', sans-serif;
            font-variation-settings: 'opsz' 64;
            font-size: 2.4em;
            font-weight: 600;
            letter-spacing: -0.02em;
            line-height: 1.05;
            color: var(--ink);
            margin: 0;
        }}
        .header p {{
            font-size: 1.05em;
            color: var(--ink-mute);
            margin-top: 12px;
            max-width: 62ch;
        }}
        .header code {{
            background: var(--surface-soft);
            border: 1px solid var(--border-soft);
            color: var(--ink);
            padding: 1px 8px;
            border-radius: 3px;
            font-size: 0.88em;
        }}

        /* Content area */
        .content {{
            padding: 0 56px 40px;
            flex: 1;
            min-width: 0;
        }}

        /* Welcome / intro card — Seqera Box style */
        .welcome-section {{
            background: var(--paper);
            border: 1px solid var(--border);
            border-radius: 3px;
            padding: 32px 36px;
            margin: 32px 0 12px;
        }}
        .welcome-section h2 {{
            font-family: 'Bricolage Grotesque', 'Inter', sans-serif;
            font-variation-settings: 'opsz' 32;
            color: var(--ink);
            font-size: 1.5em;
            font-weight: 600;
            letter-spacing: -0.01em;
            line-height: 1.2;
            margin-bottom: 14px;
        }}
        .welcome-section p {{
            color: var(--ink-mute);
            line-height: 1.65;
            margin-bottom: 12px;
        }}
        .welcome-section p:last-child {{ margin-bottom: 0; }}
        .welcome-section code {{
            background: var(--surface-soft);
            border: 1px solid var(--border-soft);
            padding: 1px 8px;
            border-radius: 3px;
            font-size: 0.88em;
            color: var(--navy);
        }}

        .section-heading {{
            font-family: 'Bricolage Grotesque', 'Inter', sans-serif;
            font-variation-settings: 'opsz' 32;
            color: var(--ink);
            margin: 48px 0 18px;
            font-size: 1.6em;
            font-weight: 600;
            letter-spacing: -0.015em;
            line-height: 1.15;
            scroll-margin-top: 24px;
            position: relative;
            padding-left: 14px;
        }}
        .section-heading::before {{
            content: "";
            position: absolute;
            left: 0;
            top: 0.3em;
            height: 0.75em;
            width: 3px;
            background: var(--teal);
            border-radius: 1px;
        }}

        /* Form fields — Seqera Box style */
        .form-group {{
            margin-bottom: 10px;
            border: 1px solid var(--border);
            border-radius: 3px;
            padding: 20px 24px;
            background: var(--paper);
            transition: border-color 0.15s ease;
        }}
        .form-group:hover {{
            border-color: rgba(0, 0, 0, 0.25);
        }}
        .form-group:focus-within {{
            border-color: var(--navy);
        }}
        .form-group.hidden {{ display: none; }}
        .form-group label.field-label {{
            display: block;
            font-family: 'JetBrains Mono', ui-monospace, monospace;
            font-weight: 500;
            color: var(--ink);
            margin-bottom: 6px;
            font-size: 0.95em;
        }}
        .form-group .description {{
            color: var(--ink-mute);
            font-size: 0.9em;
            margin-bottom: 12px;
            line-height: 1.5;
            max-width: 70ch;
        }}
        .default-badge {{
            background: var(--surface-soft);
            border: 1px solid var(--border-soft);
            color: var(--ink-mute);
            padding: 2px 8px;
            border-radius: 3px;
            font-size: 0.72em;
            margin-left: 10px;
            font-weight: 500;
            font-family: 'JetBrains Mono', ui-monospace, monospace;
            letter-spacing: 0;
            text-transform: none;
            vertical-align: middle;
        }}
        .form-group input, .form-group select, .form-group textarea {{
            width: 100%;
            padding: 9px 12px;
            border: 1px solid var(--border);
            border-radius: 4px;
            font-family: inherit;
            font-size: 0.95em;
            color: var(--ink);
            background: var(--paper);
            transition: border-color 0.15s ease, box-shadow 0.15s ease;
        }}
        .form-group input:focus, .form-group select:focus, .form-group textarea:focus {{
            outline: none;
            border-color: var(--teal-border);
            box-shadow: 0 0 0 3px rgba(49, 201, 172, 0.2);
        }}
        .form-group select {{
            appearance: none;
            -webkit-appearance: none;
            background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='12' height='12' viewBox='0 0 12 12' fill='none'%3E%3Cpath d='M3 4.5l3 3 3-3' stroke='%235c5767' stroke-width='1.5' stroke-linecap='round' stroke-linejoin='round'/%3E%3C/svg%3E");
            background-repeat: no-repeat;
            background-position: right 12px center;
            padding-right: 36px;
        }}

        .no-results {{
            text-align: center;
            padding: 60px 20px;
            color: var(--ink-mute);
            display: none;
            font-size: 0.95em;
        }}

        /* When previewing YAML, hide the editor UI */
        body.preview-mode #welcomeSection,
        body.preview-mode #sections,
        body.preview-mode #noResults {{
            display: none !important;
        }}
        body.preview-mode .section-nav a {{ pointer-events: none; opacity: 0.5; }}
        body.preview-mode .sidebar-search {{ opacity: 0.4; pointer-events: none; }}

        .btn {{
            padding: 8px 20px;
            border: 1px solid transparent;
            border-radius: 4px;
            font-family: inherit;
            font-size: 0.92em;
            cursor: pointer;
            transition: background 0.15s ease, border-color 0.15s ease, color 0.15s ease;
            font-weight: 500;
            line-height: 1.4;
        }}
        .btn-primary {{
            background: var(--navy);
            border-color: var(--navy);
            color: white;
        }}
        .btn-primary:hover {{
            background: var(--navy-2);
            border-color: var(--navy-2);
        }}
        .btn-secondary {{
            background: var(--paper);
            color: var(--ink);
            border-color: var(--border);
        }}
        .btn-secondary:hover {{
            background: var(--surface-soft);
            border-color: var(--ink-soft);
        }}
        .btn-success {{
            background: var(--teal);
            border-color: var(--teal-border);
            color: var(--navy);
        }}
        .btn-success:hover {{
            background: var(--teal-border);
            border-color: var(--teal-border-active);
            color: white;
        }}

        .yaml-output {{
            background: var(--navy);
            color: #d9d4e0;
            padding: 20px 22px;
            border-radius: 3px;
            font-family: 'JetBrains Mono', ui-monospace, SFMono-Regular, Monaco, monospace;
            font-size: 0.88em;
            line-height: 1.6;
            white-space: pre-wrap;
            max-height: 420px;
            overflow-y: auto;
            margin-top: 16px;
        }}

        .yaml-preview-section {{
            background: var(--paper);
            border: 1px solid var(--border);
            border-radius: 3px;
            padding: 28px 32px;
            margin: 32px 0;
        }}
        .yaml-preview-header {{
            margin-bottom: 18px;
        }}
        .back-btn {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 6px 14px;
            background: transparent;
            color: var(--ink-mute);
            border: 1px solid var(--border);
            border-radius: 4px;
            font-family: inherit;
            font-size: 0.88em;
            cursor: pointer;
            transition: background 0.15s ease, color 0.15s ease, border-color 0.15s ease;
        }}
        .back-btn:hover {{
            background: var(--surface-soft);
            color: var(--ink);
            border-color: var(--ink-soft);
        }}
        .back-btn span {{ font-size: 1.05em; line-height: 1; }}
        .yaml-preview-section h3 {{
            font-family: 'Bricolage Grotesque', 'Inter', sans-serif;
            font-variation-settings: 'opsz' 24;
            color: var(--ink);
            font-size: 1.2em;
            font-weight: 600;
            margin-bottom: 6px;
            letter-spacing: -0.01em;
        }}
        .yaml-preview-section p {{
            color: var(--ink-mute);
            font-size: 0.92em;
        }}
        .yaml-preview-section code {{
            background: var(--surface-soft);
            border: 1px solid var(--border-soft);
            padding: 1px 6px;
            border-radius: 3px;
            font-size: 0.85em;
        }}
        .fields-set-count {{
            color: var(--ink-soft);
            font-size: 0.82em;
            margin-top: 4px;
        }}

        @media (max-width: 900px) {{
            .sidebar {{
                position: static;
                width: 100%;
                height: auto;
            }}
            .sidebar-sections {{ max-height: none; }}
            .main {{ margin-left: 0; }}
            .header {{ padding: 32px 24px; }}
            .header h1 {{ font-size: 1.85em; }}
            .content {{ padding: 0 24px 32px; }}
            .actions {{ padding: 14px 24px; }}
            .form-group {{ padding: 16px 18px; }}
            .welcome-section {{ padding: 24px; }}
            .yaml-preview-section {{ padding: 22px 20px; }}
        }}
    </style>
</head>
<body>
    <aside class="sidebar">
        <div class="sidebar-brand">
            {MULTIQC_LOGO_SVG}
            <div class="sidebar-eyebrow">Config Wizard</div>
        </div>
        <div class="sidebar-search">
            <input type="text" id="searchInput" placeholder="Search options..." autocomplete="off" />
        </div>
        <nav class="sidebar-sections">
            <h3>Sections</h3>
            <ul class="section-nav" id="sectionNav"></ul>
        </nav>
        <div class="sidebar-actions">
            <button class="btn btn-primary" onclick="previewYaml()">Preview YAML</button>
            <button class="btn btn-success" onclick="downloadConfig()">Download Config</button>
            <button class="btn btn-secondary" onclick="resetConfig()">Reset All</button>
        </div>
        <div class="sidebar-links">
            <a href="https://seqera.io/multiqc/" target="_blank" rel="noopener">seqera.io/multiqc</a>
            <a href="https://docs.seqera.io/multiqc" target="_blank" rel="noopener">Documentation</a>
            <a href="https://github.com/MultiQC/MultiQC" target="_blank" rel="noopener">GitHub</a>
        </div>
    </aside>

    <main class="main">
        <header class="header">
            <h1>MultiQC Configuration Wizard</h1>
            <p>Build your <code>multiqc_config.yaml</code></p>
        </header>

        <div class="content" id="contentArea">
            <div class="welcome-section" id="welcomeSection">
                <h2>Welcome to the MultiQC Config Wizard</h2>
                <p>
                    This tool helps you create a custom <code>multiqc_config.yaml</code>
                    configuration file. You don't need to set every option, only
                    change what you want to customise.  Unset options will use their
                    defaults.
                </p>
                <p>
                    Browse sections in the sidebar, or use the search bar above to find
                    a specific option.  When you are done, click <strong>Preview YAML</strong>
                    or <strong>Download Config</strong> below.
                    For full documentation, visit
                    <a href="https://docs.seqera.io/multiqc" target="_blank" rel="noopener">docs.seqera.io/multiqc</a>.
                </p>
            </div>

            <div id="sections"></div>

            <div class="no-results" id="noResults">
                No configuration options match your search.
            </div>

            <div class="yaml-preview-section" id="yamlSection" style="display:none;">
                <div class="yaml-preview-header">
                    <button class="btn btn-secondary back-btn" onclick="backToEditor()">
                        <span aria-hidden="true">&larr;</span> Back to editor
                    </button>
                </div>
                <h3>Generated Configuration</h3>
                <p>Copy this content to <code>multiqc_config.yaml</code> in your project directory.</p>
                <p class="fields-set-count" id="fieldsSetCount"></p>
                <div class="yaml-output" id="yamlOutput"></div>
            </div>
        </div>
    </main>

    <script
        src="https://cdnjs.cloudflare.com/ajax/libs/js-yaml/4.1.0/js-yaml.min.js"
        integrity="sha512-CSBhVREyzHAjAFfBlIBakjoRUKp5h7VSweP0InR/pAJyptH7peuhCsqAI/snV+TwZmXZqoUklpXp6R6wMnYf5Q=="
        crossorigin="anonymous"
        referrerpolicy="no-referrer"
    ></script>
    <script>
        const configData = {config_json_escaped};

        let currentConfig = {{}};

        function initWizard() {{
            createSectionNavigation();
            createSectionContent();
            document.getElementById('searchInput').addEventListener('input', handleSearch);
        }}

        function createSectionNavigation() {{
            const nav = document.getElementById('sectionNav');
            Object.entries(configData).forEach(([sectionName, sectionProps]) => {{
                const li = document.createElement('li');
                const a = document.createElement('a');
                a.href = '#section-' + slugify(sectionName);
                a.textContent = sectionName;
                a.dataset.section = sectionName;
                const badge = document.createElement('span');
                badge.className = 'badge';
                badge.textContent = Object.keys(sectionProps).length;
                a.appendChild(badge);
                a.onclick = (e) => {{
                    e.preventDefault();
                    document.getElementById('section-' + slugify(sectionName))
                        .scrollIntoView({{ behavior: 'smooth', block: 'start' }});
                }};
                li.appendChild(a);
                nav.appendChild(li);
            }});
        }}

        function createSectionContent() {{
            const sectionsContainer = document.getElementById('sections');

            Object.entries(configData).forEach(([sectionName, sectionProps]) => {{
                const h2 = document.createElement('h2');
                h2.className = 'section-heading';
                h2.id = 'section-' + slugify(sectionName);
                h2.textContent = sectionName;
                sectionsContainer.appendChild(h2);

                Object.entries(sectionProps).forEach(([propName, propData]) => {{
                    const formGroup = createFormGroup(propName, propData);
                    formGroup.dataset.section = sectionName;
                    sectionsContainer.appendChild(formGroup);
                }});
            }});
        }}

        function createFormGroup(propName, propData) {{
            const formGroup = document.createElement('div');
            formGroup.className = 'form-group';
            formGroup.dataset.propname = propName;

            const label = document.createElement('label');
            label.className = 'field-label';
            label.textContent = propName;

            if (propData.default !== null && propData.default !== undefined) {{
                const defaultSpan = document.createElement('span');
                defaultSpan.className = 'default-badge';
                defaultSpan.textContent = 'default: ' + formatDefaultValue(propData.default);
                label.appendChild(defaultSpan);
            }}

            formGroup.appendChild(label);

            if (propData.description) {{
                const desc = document.createElement('div');
                desc.className = 'description';
                desc.textContent = propData.description;
                formGroup.appendChild(desc);
            }}

            const input = createInput(propName, propData);
            formGroup.appendChild(input);

            return formGroup;
        }}

        function createInput(propName, propData) {{
            let input;

            if (propData.enum) {{
                input = document.createElement('select');
                const emptyOption = document.createElement('option');
                emptyOption.value = '';
                emptyOption.textContent = '-- Not set --';
                input.appendChild(emptyOption);
                propData.enum.forEach(value => {{
                    const option = document.createElement('option');
                    option.value = value;
                    option.textContent = value;
                    input.appendChild(option);
                }});
                input.id = propName;
                input.onchange = () => updateConfig(propName, input);
                return input;

            }} else if (propData.type === 'boolean') {{
                input = document.createElement('select');
                input.id = propName;

                const notSet = document.createElement('option');
                notSet.value = '';
                const defaultStr = (propData.default === true) ? 'true' : (propData.default === false) ? 'false' : 'none';
                notSet.textContent = '-- Not set (default: ' + defaultStr + ') --';
                input.appendChild(notSet);

                const optTrue = document.createElement('option');
                optTrue.value = 'true';
                optTrue.textContent = 'true';
                input.appendChild(optTrue);

                const optFalse = document.createElement('option');
                optFalse.value = 'false';
                optFalse.textContent = 'false';
                input.appendChild(optFalse);

                input.onchange = () => updateConfig(propName, input);
                return input;

            }} else if (propData.type === 'array') {{
                input = document.createElement('textarea');
                input.rows = 3;
                input.placeholder = 'Enter values separated by commas or one per line';
            }} else if (propData.type === 'integer') {{
                input = document.createElement('input');
                input.type = 'number';
            }} else if (propData.type === 'number') {{
                input = document.createElement('input');
                input.type = 'number';
                input.step = 'any';
            }} else {{
                input = document.createElement('input');
                input.type = 'text';
            }}

            input.id = propName;
            input.onchange = () => updateConfig(propName, input);
            input.oninput = () => updateConfig(propName, input);
            return input;
        }}

        function updateConfig(propName, input) {{
            const propData = getCurrentPropData(propName);

            if (propData.type === 'boolean') {{
                if (input.value === '') {{
                    delete currentConfig[propName];
                }} else {{
                    currentConfig[propName] = (input.value === 'true');
                }}
            }} else if (propData.type === 'array') {{
                const value = input.value.trim();
                if (value) {{
                    currentConfig[propName] = value.split(/[,\\n]/).map(s => s.trim()).filter(s => s);
                }} else {{
                    delete currentConfig[propName];
                }}
            }} else if (propData.type === 'integer') {{
                const value = parseInt(input.value);
                if (!isNaN(value)) {{
                    currentConfig[propName] = value;
                }} else {{
                    delete currentConfig[propName];
                }}
            }} else if (propData.type === 'number') {{
                const value = parseFloat(input.value);
                if (!isNaN(value)) {{
                    currentConfig[propName] = value;
                }} else {{
                    delete currentConfig[propName];
                }}
            }} else {{
                const value = input.value.trim();
                if (value) {{
                    currentConfig[propName] = value;
                }} else {{
                    delete currentConfig[propName];
                }}
            }}
        }}

        function getCurrentPropData(propName) {{
            for (const sectionProps of Object.values(configData)) {{
                if (sectionProps[propName]) return sectionProps[propName];
            }}
            return {{}};
        }}

        function handleSearch() {{
            const query = document.getElementById('searchInput').value.toLowerCase().trim();
            const formGroups = document.querySelectorAll('.form-group');
            const headings = document.querySelectorAll('.section-heading');
            let anyVisible = false;

            formGroups.forEach(fg => {{
                const name = (fg.dataset.propname || '').toLowerCase();
                const desc = (fg.querySelector('.description') || {{}}).textContent || '';
                if (!query || name.includes(query) || desc.toLowerCase().includes(query)) {{
                    fg.classList.remove('hidden');
                    anyVisible = true;
                }} else {{
                    fg.classList.add('hidden');
                }}
            }});

            headings.forEach(h => {{
                const sectionSlug = h.id;
                let sibling = h.nextElementSibling;
                let hasVisible = false;
                while (sibling && !sibling.classList.contains('section-heading')) {{
                    if (sibling.classList.contains('form-group') && !sibling.classList.contains('hidden')) {{
                        hasVisible = true;
                    }}
                    sibling = sibling.nextElementSibling;
                }}
                h.style.display = (!query || hasVisible) ? '' : 'none';
            }});

            document.getElementById('noResults').style.display = (query && !anyVisible) ? 'block' : 'none';
            document.getElementById('welcomeSection').style.display = query ? 'none' : '';
        }}

        function formatDefaultValue(value) {{
            if (Array.isArray(value)) {{
                if (value.length <= 3) {{
                    return JSON.stringify(value);
                }}
                return '[' + value.length + ' items]';
            }} else if (typeof value === 'object' && value !== null) {{
                const keys = Object.keys(value);
                if (keys.length <= 2) {{
                    return JSON.stringify(value);
                }}
                return '{{' + keys.length + ' keys}}';
            }} else if (typeof value === 'string') {{
                return '"' + value + '"';
            }}
            return String(value);
        }}

        function generateYaml() {{
            const output = document.getElementById('yamlOutput');
            const count = Object.keys(currentConfig).length;
            document.getElementById('fieldsSetCount').textContent = count + ' option' + (count !== 1 ? 's' : '') + ' configured';

            if (count === 0) {{
                output.textContent = '# No configuration set yet\\n# Change some values above, then come back here.';
                return;
            }}

            try {{
                output.textContent = jsyaml.dump(currentConfig, {{
                    indent: 2,
                    lineWidth: 80,
                    noRefs: true,
                    sortKeys: true
                }});
            }} catch (error) {{
                output.textContent = '# Error generating YAML: ' + error.message;
            }}
        }}

        function previewYaml() {{
            generateYaml();
            document.body.classList.add('preview-mode');
            document.getElementById('yamlSection').style.display = 'block';
            window.scrollTo({{ top: 0, behavior: 'smooth' }});
        }}

        function backToEditor() {{
            document.body.classList.remove('preview-mode');
            document.getElementById('yamlSection').style.display = 'none';
            window.scrollTo({{ top: 0, behavior: 'smooth' }});
        }}

        function downloadConfig() {{
            if (Object.keys(currentConfig).length === 0) {{
                alert('No configuration to download. Please set some options first.');
                return;
            }}
            generateYaml();
            const yaml = document.getElementById('yamlOutput').textContent;
            const blob = new Blob([yaml], {{ type: 'text/yaml' }});
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'multiqc_config.yaml';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }}

        function resetConfig() {{
            if (!confirm('Reset all configuration? This cannot be undone.')) return;
            currentConfig = {{}};
            document.querySelectorAll('.form-group input, .form-group select, .form-group textarea').forEach(input => {{
                if (input.tagName === 'SELECT') {{
                    input.selectedIndex = 0;
                }} else {{
                    input.value = '';
                }}
            }});
            if (document.body.classList.contains('preview-mode')) {{
                generateYaml();
            }}
        }}

        function slugify(text) {{
            return text.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/(^-|-$)/g, '');
        }}

        document.addEventListener('DOMContentLoaded', function() {{
            if (typeof jsyaml === 'undefined') {{
                document.getElementById('yamlOutput').textContent =
                    '# Error: js-yaml library failed to load.\\n# Please check your internet connection and refresh.';
            }} else {{
                initWizard();
            }}
        }});
    </script>
</body>
</html>"""


if __name__ == "__main__":
    html_output = generate_config_wizard()

    output_path = Path(__file__).parent.parent / "docs" / "multiqc_config_wizard.html"

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(html_output)
    except OSError as e:
        print(f"Error writing to {output_path}: {e}")
        sys.exit(1)

    print(f"MultiQC Configuration Wizard generated at {output_path}")
    print("Open the HTML file in your web browser to use the wizard.")
