# MultiQC Dev Tools Plugin Marketplace

A curated collection of Claude Code plugins for bioinformatics tool developers, focused on MultiQC development workflows and similar projects.

## Overview

This marketplace provides specialized skills that help Claude Code assist with common bioinformatics development tasks:

- **Module Generation**: Automatically create complete MultiQC modules from specifications
- **Log Analysis**: Understand and parse unfamiliar bioinformatics tool log files
- **Template Development**: Build custom MultiQC report templates

These plugins are designed for the bioinformatics community but can be adapted for similar data analysis and reporting tools.

## Installation

### Add the Marketplace

In your Claude Code session, add this marketplace:

```bash
/plugin marketplace add https://github.com/MultiQC/MultiQC/.claude-plugin/marketplace.json
```

Or if working locally with the MultiQC repository:

```bash
/plugin marketplace add /path/to/MultiQC/.claude-plugin/marketplace.json
```

### Install Individual Plugins

After adding the marketplace, install the plugins you need:

```bash
# Install all plugins
/plugin install creating-multiqc-modules
/plugin install analyzing-log-files
/plugin install creating-multiqc-template

# Or list available plugins first
/plugin marketplace list
```

## Available Plugins

### 1. Creating MultiQC Modules

**Plugin**: `creating-multiqc-modules`

**Purpose**: Generate complete, standards-compliant MultiQC module implementations from GitHub issue specifications or example log files.

**What it does**:

- Downloads and analyzes example log files from issues
- Identifies data structures and optimal parsing strategies
- Generates production-ready module code following MultiQC best practices
- Creates comprehensive unit tests
- Updates configuration files (search_patterns.yaml, pyproject.toml)
- Prepares pull requests with proper formatting

**When to use**:

- Creating new MultiQC modules from scratch
- Converting example output files into working code
- Implementing module requests from GitHub issues
- Need to follow MultiQC coding standards automatically

**Example usage**:

```
You: "Create a module for the FastQC tool based on issue #1234"
Claude: [Automatically uses the module creation skill to generate complete implementation]
```

**Included resources**:

- Production-ready module template
- Examples of simple and complex modules
- Supporting documentation for parsing strategies

---

### 2. Analyzing Log Files

**Plugin**: `analyzing-log-files`

**Purpose**: Systematically examine unfamiliar bioinformatics tool output files to understand their structure and design optimal parsing approaches.

**What it does**:

- Identifies file format (JSON, TSV, key-value pairs, custom formats)
- Detects section structure and delimiter patterns
- Locates sample name extraction points
- Identifies key metrics and their data types
- Suggests appropriate parsing strategies (regex, JSON parsing, etc.)
- Recommends visualization types for the data
- Plans which metrics should appear in general statistics

**When to use**:

- Encountering new or unfamiliar log file formats
- Designing parsing logic for MultiQC modules
- Analyzing example files from module requests
- Troubleshooting parsing problems in existing modules

**Example usage**:

```
You: "I need to parse this bcftools stats output but I'm not sure of the best approach"
Claude: [Analyzes the file structure and recommends parsing strategy with code examples]
```

**Included resources**:

- Common pattern reference guide
- Examples of different file format types
- Parsing strategy decision trees

---

### 3. Creating MultiQC Templates

**Plugin**: `creating-multiqc-template`

**Purpose**: Guide developers through creating custom MultiQC report templates for branded or specialized reporting needs.

**What it does**:

- Explains MultiQC template architecture and structure
- Guides Python package setup for templates
- Helps customize Jinja2 template files
- Manages static assets (CSS, JavaScript, images)
- Implements custom plotting functions
- Handles child template inheritance patterns

**When to use**:

- Building custom report themes for your organization
- Creating specialized report layouts
- Customizing report branding and styling
- Adding custom interactive visualizations
- Developing reusable template components

**Example usage**:

```
You: "I need to create a custom MultiQC template with my company's branding"
Claude: [Guides through template creation with proper structure and examples]
```

**Included resources**:

- Template structure guide
- Jinja2 customization examples
- Asset management best practices
- Child template patterns

## Target Audience

These plugins are designed for:

- **Bioinformatics tool developers** creating analysis software
- **Core facility developers** building custom reporting solutions
- **MultiQC contributors** developing new modules
- **Computational biologists** integrating tools into MultiQC
- **DevOps teams** customizing MultiQC for organizational needs

## Requirements

- Claude Code CLI or IDE integration
- Basic familiarity with:
  - Python programming
  - Bioinformatics file formats
  - MultiQC architecture (for module and template plugins)

## Why Use These Plugins?

### Automated Best Practices

These plugins encode years of MultiQC development experience, ensuring:

- Consistent code style and structure
- Proper error handling patterns
- Comprehensive test coverage
- Standards-compliant implementations

### Faster Development

- **Module creation**: Hours → Minutes
- **Log parsing**: Trial and error → Systematic approach
- **Template development**: Days → Hours

### Knowledge Transfer

Embedded documentation and examples help developers learn:

- MultiQC architecture patterns
- Python best practices for data analysis
- Effective visualization strategies
- Testing and quality assurance approaches

## Examples

### Complete Workflow: New Module Request

```
1. You receive a GitHub issue requesting a new module

You: "Analyze issue #5678 for creating a module for tool XYZ"

2. Claude uses analyzing-log-files to understand the format
   - Identifies TSV structure with headers
   - Locates sample names in first column
   - Finds key metrics for visualization

3. Claude uses creating-multiqc-modules to generate code
   - Creates complete module implementation
   - Writes unit tests
   - Updates configuration files
   - Prepares pull request

Result: Production-ready module in minutes
```

### Customization Workflow

```
You: "I need to create a branded report template for our core facility"

Claude uses creating-multiqc-template to:
1. Set up proper Python package structure
2. Create base template files
3. Add custom CSS with your branding
4. Implement specialized plot types
5. Set up inheritance from default template

Result: Custom template ready to deploy
```

## Support and Contributions

### Getting Help

If you encounter issues:

1. Check the included documentation in each plugin
2. Review examples in `.claude/skills/*/examples/`
3. Consult the MultiQC documentation
4. Open an issue on the MultiQC GitHub repository

### Contributing

To improve these plugins:

1. Fork the MultiQC repository
2. Make changes to skills in `.claude/skills/`
3. Test thoroughly with multiple scenarios
4. Submit a pull request with clear descriptions

### Feedback

We welcome feedback on:

- Plugin effectiveness and accuracy
- Missing features or capabilities
- Documentation clarity
- Use cases we haven't considered

Open an issue or discussion in the MultiQC repository.

## Technical Details

### Plugin Structure

Each plugin is a Claude Code skill with:

- `SKILL.md`: Main skill definition with YAML frontmatter
- Supporting documentation files
- Code templates and examples
- Reference materials

### Source Control

Plugins reference the `.claude/skills/` directory in the MultiQC repository. Updates to skills automatically propagate to plugin users.

### Versioning

- Plugins follow semantic versioning (MAJOR.MINOR.PATCH)
- Breaking changes increment MAJOR version
- New features increment MINOR version
- Bug fixes increment PATCH version

### Privacy and Security

These plugins:

- Run entirely locally in Claude Code
- Do not send data to external services
- Only access files you explicitly provide
- Follow Claude Code's security model

## License

These plugins are distributed under the MIT License, consistent with the MultiQC project.

## Acknowledgments

Developed by the MultiQC team with contributions from the bioinformatics community.

Special thanks to:

- MultiQC contributors and maintainers
- Claude Code team at Anthropic
- Bioinformatics developers who provided feedback

---

**Last Updated**: 2025-10-27
**Version**: 1.0.0
**Maintained by**: MultiQC Development Team
**Repository**: https://github.com/MultiQC/MultiQC
