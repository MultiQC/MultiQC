#!/usr/bin/env python3
"""
Script to generate a web-based configuration wizard for MultiQC.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml

# Add parent directory to path so we can import multiqc
sys.path.insert(0, str(Path(__file__).parent.parent))

from multiqc.utils.config_schema import MultiQCConfig


def generate_config_wizard():
    """Generate a web-based configuration wizard for MultiQC."""

    # Get JSON schema and default values
    schema = MultiQCConfig.model_json_schema()
    properties = schema.get("properties", {})

    # Load default values from config_defaults.yaml
    config_defaults_path = Path(__file__).parent.parent / "multiqc" / "config_defaults.yaml"
    with open(config_defaults_path, "r") as f:
        config_defaults = yaml.safe_load(f)

    # Group properties into logical sections
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
            "sample_names_ignore",
            "sample_names_ignore_re",
            "sample_names_only_include",
            "sample_names_only_include_re",
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
        "File Discovery": ["require_logs", "ignore_symlinks", "ignore_images", "fn_ignore_dirs", "fn_ignore_paths"],
    }

    # Build configuration data for JavaScript
    config_data = {}
    for section_name, section_props in sections.items():
        config_data[section_name] = {}
        for prop_name in section_props:
            if prop_name in properties:
                prop = properties[prop_name]

                # Get type information
                prop_type = prop.get("type", "string")
                if "anyOf" in prop:
                    # Handle Optional types
                    type_options = [t.get("type") for t in prop.get("anyOf", []) if "type" in t]
                    if type_options:
                        prop_type = type_options[0]

                # Handle enum/literal types
                enum_values = None
                if "enum" in prop:
                    enum_values = prop["enum"]
                elif "$defs" in schema and "AiProviderLiteral" in schema["$defs"]:
                    if prop_name == "ai_provider":
                        enum_values = ["seqera", "openai", "anthropic", "aws_bedrock", "custom"]

                config_data[section_name][prop_name] = {
                    "type": prop_type,
                    "description": prop.get("description", ""),
                    "default": config_defaults.get(prop_name),
                    "enum": enum_values,
                }

    # Generate HTML
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MultiQC Configuration Wizard</title>
    <style>
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 12px;
            box-shadow: 0 20px 40px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #2c3e50 0%, #3498db 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        
        .main-content {{
            display: flex;
            min-height: 70vh;
        }}
        
        .sidebar {{
            width: 300px;
            background: #f8f9fa;
            border-right: 1px solid #e9ecef;
            padding: 20px;
        }}
        
        .sidebar h3 {{
            color: #2c3e50;
            margin-bottom: 15px;
            font-size: 1.2em;
        }}
        
        .section-nav {{
            list-style: none;
        }}
        
        .section-nav li {{
            margin-bottom: 8px;
        }}
        
        .section-nav a {{
            color: #6c757d;
            text-decoration: none;
            padding: 8px 12px;
            border-radius: 6px;
            display: block;
            transition: all 0.3s ease;
        }}
        
        .section-nav a:hover, .section-nav a.active {{
            background: #3498db;
            color: white;
        }}
        
        .content {{
            flex: 1;
            padding: 30px;
            overflow-y: auto;
        }}
        
        .section {{
            display: none;
            animation: fadeIn 0.3s ease;
        }}
        
        .section.active {{
            display: block;
        }}
        
        @keyframes fadeIn {{
            from {{ opacity: 0; transform: translateY(10px); }}
            to {{ opacity: 1; transform: translateY(0); }}
        }}
        
        .section h2 {{
            color: #2c3e50;
            margin-bottom: 20px;
            font-size: 1.8em;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
        }}
        
        .form-group {{
            margin-bottom: 25px;
            border: 1px solid #e9ecef;
            border-radius: 8px;
            padding: 20px;
            background: #fafbfc;
            transition: all 0.3s ease;
        }}
        
        .form-group:hover {{
            border-color: #3498db;
            box-shadow: 0 2px 8px rgba(52, 152, 219, 0.1);
        }}
        
        .form-group label {{
            display: block;
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 8px;
            font-size: 1.1em;
        }}
        
        .form-group .description {{
            color: #6c757d;
            font-size: 0.9em;
            margin-bottom: 12px;
            line-height: 1.4;
        }}
        
        .form-group input, .form-group select, .form-group textarea {{
            width: 100%;
            padding: 12px;
            border: 1px solid #ddd;
            border-radius: 6px;
            font-size: 1em;
            transition: border-color 0.3s ease;
        }}
        
        .form-group input:focus, .form-group select:focus, .form-group textarea:focus {{
            outline: none;
            border-color: #3498db;
            box-shadow: 0 0 0 3px rgba(52, 152, 219, 0.1);
        }}
        
        .form-group input[type="checkbox"] {{
            width: auto;
            margin-right: 8px;
        }}
        
        .checkbox-wrapper {{
            display: flex;
            align-items: center;
        }}
        
        .default-value {{
            background: #e8f4fd;
            color: #2c5282;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            margin-left: 8px;
        }}
        
        .actions {{
            position: sticky;
            bottom: 0;
            background: white;
            border-top: 1px solid #e9ecef;
            padding: 20px 30px;
            display: flex;
            gap: 15px;
            justify-content: space-between;
            align-items: center;
        }}
        
        .btn {{
            padding: 12px 24px;
            border: none;
            border-radius: 6px;
            font-size: 1em;
            cursor: pointer;
            text-decoration: none;
            display: inline-block;
            transition: all 0.3s ease;
            font-weight: 600;
        }}
        
        .btn-primary {{
            background: #3498db;
            color: white;
        }}
        
        .btn-primary:hover {{
            background: #2980b9;
            transform: translateY(-1px);
        }}
        
        .btn-secondary {{
            background: #6c757d;
            color: white;
        }}
        
        .btn-secondary:hover {{
            background: #5a6268;
        }}
        
        .btn-success {{
            background: #27ae60;
            color: white;
        }}
        
        .btn-success:hover {{
            background: #229954;
        }}
        
        .progress {{
            height: 4px;
            background: #e9ecef;
            border-radius: 2px;
            overflow: hidden;
        }}
        
        .progress-bar {{
            height: 100%;
            background: #3498db;
            transition: width 0.3s ease;
        }}
        
        .yaml-output {{
            background: #2c3e50;
            color: #ecf0f1;
            padding: 20px;
            border-radius: 8px;
            font-family: 'Monaco', 'Courier New', monospace;
            font-size: 0.9em;
            white-space: pre-wrap;
            max-height: 400px;
            overflow-y: auto;
            margin-top: 20px;
        }}
        
        @media (max-width: 768px) {{
            .main-content {{
                flex-direction: column;
            }}
            
            .sidebar {{
                width: 100%;
                border-right: none;
                border-bottom: 1px solid #e9ecef;
            }}
            
            .section-nav {{
                display: flex;
                flex-wrap: wrap;
                gap: 5px;
            }}
            
            .section-nav li {{
                margin-bottom: 0;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 MultiQC Configuration Wizard</h1>
            <p>Create your perfect MultiQC configuration with guided assistance</p>
        </div>
        
        <div class="main-content">
            <div class="sidebar">
                <h3>Configuration Sections</h3>
                <ul class="section-nav" id="sectionNav">
                    <!-- Will be populated by JavaScript -->
                </ul>
                
                <div class="progress" style="margin-top: 20px;">
                    <div class="progress-bar" id="progressBar" style="width: 0%"></div>
                </div>
                <p style="margin-top: 10px; color: #6c757d; font-size: 0.9em;">
                    <span id="progressText">0% complete</span>
                </p>
            </div>
            
            <div class="content">
                <div id="sections">
                    <!-- Will be populated by JavaScript -->
                </div>
                
                <div id="yamlSection" class="section">
                    <h2>🎉 Generated Configuration</h2>
                    <p>Here's your custom MultiQC configuration file. Copy this content to <code>multiqc_config.yaml</code> in your project directory.</p>
                    <div class="yaml-output" id="yamlOutput"></div>
                </div>
            </div>
        </div>
        
        <div class="actions">
            <div>
                <button class="btn btn-secondary" onclick="resetConfig()">🔄 Reset All</button>
                <button class="btn btn-secondary" onclick="loadDefaults()">⚙️ Load Defaults</button>
            </div>
            <div>
                <button class="btn btn-primary" onclick="previewYaml()">👁️ Preview YAML</button>
                <button class="btn btn-success" onclick="downloadConfig()">📥 Download Config</button>
            </div>
        </div>
    </div>

    <script src="https://cdnjs.cloudflare.com/ajax/libs/js-yaml/4.1.0/js-yaml.min.js"></script>
    <script>
        // Configuration data from Python
        const configData = {json.dumps(config_data, indent=8)};
        
        let currentConfig = {{}};
        let currentSection = null;
        let yamlNavLink = null; // Store reference to YAML nav link
        
        // Initialize the wizard
        function initWizard() {{
            createSectionNavigation();
            createSectionContent();
            showSection(Object.keys(configData)[0]);
        }}
        
        function createSectionNavigation() {{
            const nav = document.getElementById('sectionNav');
            Object.keys(configData).forEach((sectionName, index) => {{
                const li = document.createElement('li');
                const a = document.createElement('a');
                a.href = '#';
                a.textContent = sectionName;
                a.dataset.section = sectionName; // Add data attribute for easier selection
                a.onclick = (e) => {{
                    e.preventDefault();
                    showSection(sectionName);
                }};
                li.appendChild(a);
                nav.appendChild(li);
            }});
            
            // Add YAML preview section
            const yamlLi = document.createElement('li');
            const yamlA = document.createElement('a');
            yamlA.href = '#';
            yamlA.textContent = '📝 Preview YAML';
            yamlA.dataset.section = 'yaml'; // Add data attribute
            yamlA.onclick = (e) => {{
                e.preventDefault();
                showSection('yaml');
            }};
            yamlLi.appendChild(yamlA);
            nav.appendChild(yamlLi);
            yamlNavLink = yamlA; // Store reference
        }}
        
        function createSectionContent() {{
            const sectionsContainer = document.getElementById('sections');
            
            Object.entries(configData).forEach(([sectionName, sectionProps]) => {{
                const sectionDiv = document.createElement('div');
                sectionDiv.className = 'section';
                sectionDiv.id = sectionName.replace(/\\s+/g, '-').toLowerCase();
                
                const h2 = document.createElement('h2');
                h2.textContent = sectionName;
                sectionDiv.appendChild(h2);
                
                Object.entries(sectionProps).forEach(([propName, propData]) => {{
                    const formGroup = createFormGroup(propName, propData);
                    sectionDiv.appendChild(formGroup);
                }});
                
                sectionsContainer.appendChild(sectionDiv);
            }});
        }}
        
        function createFormGroup(propName, propData) {{
            const formGroup = document.createElement('div');
            formGroup.className = 'form-group';
            
            const label = document.createElement('label');
            label.textContent = propName;
            
            if (propData.default !== null && propData.default !== undefined) {{
                const defaultSpan = document.createElement('span');
                defaultSpan.className = 'default-value';
                defaultSpan.textContent = `default: ${{formatDefaultValue(propData.default)}}`;
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
                
                // Add empty option
                const emptyOption = document.createElement('option');
                emptyOption.value = '';
                emptyOption.textContent = '-- Select --';
                input.appendChild(emptyOption);
                
                propData.enum.forEach(value => {{
                    const option = document.createElement('option');
                    option.value = value;
                    option.textContent = value;
                    input.appendChild(option);
                }});
            }} else if (propData.type === 'boolean') {{
                const wrapper = document.createElement('div');
                wrapper.className = 'checkbox-wrapper';
                
                input = document.createElement('input');
                input.type = 'checkbox';
                input.id = propName; // Set ID on the actual input, not the wrapper
                
                const label = document.createElement('label');
                label.textContent = 'Enabled';
                label.style.fontWeight = 'normal';
                label.style.marginBottom = '0';
                label.style.cursor = 'pointer';
                label.onclick = () => {{
                    input.checked = !input.checked;
                    updateConfig(propName, input);
                }};
                
                input.onchange = () => updateConfig(propName, input);
                input.oninput = () => updateConfig(propName, input);
                
                wrapper.appendChild(input);
                wrapper.appendChild(label);
                return wrapper;
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
            
            if (propData.type !== 'boolean') {{
                input.id = propName;
                input.onchange = () => updateConfig(propName, input);
                input.oninput = () => updateConfig(propName, input);
            }}
            
            return input;
        }}
        
        function updateConfig(propName, input) {{
            const propData = getCurrentPropData(propName);
            
            if (propData.type === 'boolean') {{
                // Always set boolean values explicitly (true or false)
                currentConfig[propName] = input.checked;
            }} else if (propData.type === 'array') {{
                const value = input.value.trim();
                if (value) {{
                    // Split by comma or newline and clean up
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
            
            updateProgress();
        }}
        
        function getCurrentPropData(propName) {{
            for (const sectionProps of Object.values(configData)) {{
                if (sectionProps[propName]) {{
                    return sectionProps[propName];
                }}
            }}
            return {{}};
        }}
        
        function showSection(sectionName) {{
            // Update navigation - use data attributes for reliable selection
            document.querySelectorAll('.section-nav a').forEach(a => a.classList.remove('active'));
            
            if (sectionName === 'yaml') {{
                if (yamlNavLink) {{
                    yamlNavLink.classList.add('active');
                }}
                document.querySelectorAll('.section').forEach(s => s.classList.remove('active'));
                document.getElementById('yamlSection').classList.add('active');
                generateYaml();
            }} else {{
                const targetNav = document.querySelector(`.section-nav a[data-section="${{sectionName}}"]`);
                if (targetNav) {{
                    targetNav.classList.add('active');
                }}
                
                document.querySelectorAll('.section').forEach(s => s.classList.remove('active'));
                const targetSection = document.getElementById(sectionName.replace(/\\s+/g, '-').toLowerCase());
                if (targetSection) {{
                    targetSection.classList.add('active');
                }}
            }}
            
            currentSection = sectionName;
        }}
        
        function updateProgress() {{
            const totalFields = Object.values(configData).reduce((sum, section) => 
                sum + Object.keys(section).length, 0
            );
            const filledFields = Object.keys(currentConfig).length;
            const percentage = Math.round((filledFields / totalFields) * 100);
            
            document.getElementById('progressBar').style.width = percentage + '%';
            document.getElementById('progressText').textContent = `${{percentage}}% complete (${{filledFields}}/${{totalFields}} fields)`;
        }}
        
        function formatDefaultValue(value) {{
            if (Array.isArray(value)) {{
                return `[${{value.length}} items]`;
            }} else if (typeof value === 'object' && value !== null) {{
                return `{{${{Object.keys(value).length}} keys}}`;
            }} else if (typeof value === 'string') {{
                return `"${{value}}"`;
            }}
            return String(value);
        }}
        
        function generateYaml() {{
            try {{
                if (Object.keys(currentConfig).length === 0) {{
                    document.getElementById('yamlOutput').textContent = '# No configuration set yet\\n# Fill out some fields in the other sections and come back here to see your YAML config!';
                    return;
                }}
                
                const yaml = jsyaml.dump(currentConfig, {{
                    indent: 2,
                    lineWidth: 80,
                    noRefs: true,
                    sortKeys: true
                }});
                document.getElementById('yamlOutput').textContent = yaml;
            }} catch (error) {{
                console.error('Error generating YAML:', error);
                document.getElementById('yamlOutput').textContent = '# Error generating YAML: ' + error.message;
            }}
        }}
        
        function previewYaml() {{
            showSection('yaml');
        }}
        
        function downloadConfig() {{
            if (Object.keys(currentConfig).length === 0) {{
                alert('No configuration to download! Please fill out some fields first.');
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
            if (confirm('Are you sure you want to reset all configuration? This cannot be undone.')) {{
                currentConfig = {{}};
                document.querySelectorAll('input, select, textarea').forEach(input => {{
                    if (input.type === 'checkbox') {{
                        input.checked = false;
                    }} else {{
                        input.value = '';
                    }}
                }});
                updateProgress();
                if (currentSection === 'yaml') {{
                    generateYaml(); // Refresh YAML if currently viewing it
                }}
            }}
        }}
        
        function loadDefaults() {{
            if (confirm('Load default values for all fields? This will overwrite your current settings.')) {{
                currentConfig = {{}};
                Object.entries(configData).forEach(([sectionName, sectionProps]) => {{
                    Object.entries(sectionProps).forEach(([propName, propData]) => {{
                        if (propData.default !== null && propData.default !== undefined) {{
                            currentConfig[propName] = propData.default;
                            
                            const input = document.getElementById(propName);
                            if (input) {{
                                if (propData.type === 'boolean') {{
                                    input.checked = propData.default;
                                }} else if (propData.type === 'array') {{
                                    input.value = Array.isArray(propData.default) ? 
                                        propData.default.join(', ') : String(propData.default);
                                }} else {{
                                    input.value = propData.default;
                                }}
                            }}
                        }}
                    }});
                }});
                updateProgress();
                if (currentSection === 'yaml') {{
                    generateYaml(); // Refresh YAML if currently viewing it
                }}
            }}
        }}
        
        // Check if js-yaml loaded properly
        document.addEventListener('DOMContentLoaded', function() {{
            if (typeof jsyaml === 'undefined') {{
                console.error('js-yaml library failed to load');
                document.getElementById('yamlOutput').textContent = '# Error: js-yaml library failed to load\\n# Please check your internet connection and refresh the page.';
            }} else {{
                initWizard();
            }}
        }});
    </script>
</body>
</html>"""

    return html_content


if __name__ == "__main__":
    # Generate the wizard HTML
    html = generate_config_wizard()

    # Output file
    output_path = Path(__file__).parent.parent / "docs" / "multiqc_config_wizard.html"

    # Write to file
    with open(output_path, "w") as f:
        f.write(html)

    print(f"MultiQC Configuration Wizard generated at {output_path}")
    print("Open the HTML file in your web browser to use the wizard.")
