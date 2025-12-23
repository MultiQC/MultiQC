# MultiQC Configuration Wizard

The MultiQC Configuration Wizard is a web-based tool that helps you create custom MultiQC configuration files with guided assistance and autosuggestions.

## Features

✨ **Interactive Interface** - Easy-to-use web interface with form sections
🔍 **Auto-suggestions** - Dropdown menus for enumerated values (like AI providers)  
📖 **Built-in Documentation** - Each field shows its description and default value
📊 **Progress Tracking** - Visual progress bar showing completion status
💾 **YAML Generation** - Real-time YAML configuration file generation
📥 **Download Support** - Download your configuration as a YAML file
🔄 **Reset & Defaults** - Reset all settings or load default values
📱 **Mobile Responsive** - Works on desktop, tablet, and mobile devices

## How to Use

1. **Open the Wizard**

   ```bash
   # Generate the wizard (if not already done)
   python scripts/generate_config_wizard.py

   # Open the HTML file in your browser
   open docs/multiqc_config_wizard.html
   ```

2. **Configure MultiQC**
   - Navigate through different configuration sections using the sidebar
   - Fill in the fields you want to customize
   - See descriptions and default values for each option
   - Track your progress with the progress bar

3. **Generate Configuration**
   - Click "Preview YAML" to see your configuration
   - Click "Download Config" to save the `multiqc_config.yaml` file

4. **Use Your Configuration**
   - Place the downloaded `multiqc_config.yaml` file in your project directory
   - Run MultiQC normally - it will automatically use your configuration

## Configuration Sections

The wizard organizes MultiQC options into logical sections:

- **Report Appearance** - Title, subtitle, logos, styling
- **Output Options** - File names, formats, export settings
- **AI Summary** - AI-powered report summaries
- **Plot Settings** - Chart appearance and behavior
- **Table Settings** - Data table formatting
- **Sample Names** - Sample name cleaning and formatting
- **Performance & Debugging** - Logging and performance options
- **File Discovery** - File search and filtering options

## Tips

- **Hover over fields** to see additional styling
- **Use the progress bar** to track how much you've configured
- **Load defaults** to see typical MultiQC settings
- **Preview YAML** frequently to see your configuration taking shape
- **Reset all** if you want to start over

## Example Workflow

1. Set your report title and subtitle in "Report Appearance"
2. Configure output file names in "Output Options"
3. Enable AI summaries in "AI Summary" if desired
4. Adjust plot settings for your preferred visualization style
5. Configure sample name cleaning if needed
6. Preview and download your configuration

## Technical Details

The wizard is a self-contained HTML file that:

- Uses vanilla JavaScript (no external dependencies except js-yaml)
- Loads configuration schema from MultiQC's Pydantic models
- Includes default values from `config_defaults.yaml`
- Generates valid YAML configuration files
- Works offline once loaded

## Regenerating the Wizard

To update the wizard with new configuration options:

```bash
python scripts/generate_config_wizard.py
```

This will regenerate the HTML file with the latest configuration schema and default values.
