---
name: Creating MultiQC Template
description: Guide for developing custom MultiQC report templates to customize report appearance, branding, and visualization. Covers Python package structure, Jinja2 templates, asset management, custom plotting functions, and child template inheritance. Keywords: MultiQC, template development, report customization, Jinja2, Python submodule, custom visualization, branding, HTML reports
argument-hint: <template-name>
---

# Creating MultiQC Template

Comprehensive guide for developing custom MultiQC report templates to customize appearance, add institutional branding, or create specialized visualizations.

## Table of Contents

1. [Template Architecture Overview](#template-architecture-overview)
2. [Project Setup and Structure](#project-setup-and-structure)
3. [Template Development](#template-development)
4. [Asset Management](#asset-management)
5. [Custom Plotting](#custom-plotting)
6. [Development Workflow](#development-workflow)
7. [Testing and Validation](#testing-and-validation)
8. [Distribution](#distribution)

---

## Template Architecture Overview

### Template Types

**Standard Template:**

- Completely standalone implementation
- Full control over all aspects
- Requires implementing all functionality
- Example: `default`, `geo`, `sections` templates

**Child Template (Recommended):**

- Inherits from existing template (usually `default`)
- Override only specific components
- Leverage 1448+ lines of parent JavaScript
- Perfect for branding, minor customizations
- Example: Custom header with institutional logo

### When to Use Each Type

**Use Child Template when:**

- Adding institutional branding/logos
- Customizing colors and styling
- Modifying headers/footers
- Tweaking existing layouts
- 90% of use cases

**Use Standard Template when:**

- Completely different report structure needed
- Novel visualization approaches
- Integration with external frameworks
- Publishing as general-purpose template

---

## Project Setup and Structure

### Minimal Template Structure

```
multiqc_<template-name>/                 # Template package directory
├── setup.py                             # Package configuration
├── multiqc_<template-name>/             # Template module
│   ├── __init__.py                      # Template registration
│   ├── <template-name>.py               # Optional: custom functions
│   └── templates/                       # Template files directory
│       ├── base.html                    # Main Jinja2 template
│       ├── header.html                  # Custom header (optional)
│       ├── footer.html                  # Custom footer (optional)
│       └── assets/                      # Static resources
│           ├── css/
│           │   └── custom.css          # Custom styles
│           ├── js/
│           │   └── custom.js           # Custom JavaScript
│           └── img/
│               └── logo.png            # Images/logos
├── README.md                            # Documentation
└── LICENSE                              # License file
```

### Setup Script (setup.py)

**Standard template:**

```python
from setuptools import setup, find_packages

setup(
    name='multiqc_<template-name>',
    version='1.0.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='Custom MultiQC template for <purpose>',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/multiqc_<template-name>',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'multiqc>=1.14'
    ],
    entry_points={
        'multiqc.templates.v1': [
            '<template-name> = multiqc_<template-name>'
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    python_requires='>=3.7',
)
```

**Key components:**

- `entry_points`: Registers template with MultiQC
- `include_package_data=True`: Bundles template files
- `multiqc.templates.v1`: Entry point group name

### Template Initialization (**init**.py)

**Child template example:**

```python
#!/usr/bin/env python
"""
MultiQC template for <organization> reports.
Inherits from default template with custom branding.
"""

import os

# Inherit from default template
template_parent = 'default'

# Path to template directory
template_dir = os.path.dirname(__file__)

# Main Jinja2 template file
base_fn = 'base.html'

# Output subdirectory (optional)
# output_subdir = '<template-name>_report'

# Files to copy to output directory
copy_files = [
    'assets',  # Copy entire assets directory
]

# Configuration overrides (optional)
# config.plots_force_flat = True  # Force static images instead of interactive
```

**Standalone template example:**

```python
#!/usr/bin/env python
"""
Custom MultiQC template with novel visualization approach.
"""

import os

# NO template_parent - this is a standalone template

template_dir = os.path.dirname(__file__)
base_fn = 'base.html'

# Custom plotting functions (optional)
def bargraph(plotdata, pconfig):
    """
    Custom bar graph implementation.
    Return HTML string or None to use default.
    """
    # Your custom visualization logic here
    return None  # Return None to use default implementation

def linegraph(plotdata, pconfig):
    """
    Custom line graph implementation.
    """
    return None
```

**Key variables:**

- `template_parent`: Inherit from another template (omit for standalone)
- `template_dir`: Path to template files (always `os.path.dirname(__file__)`)
- `base_fn`: Main Jinja2 template filename
- `output_subdir`: Subdirectory for output files (optional)
- `copy_files`: List of files/directories to bundle with report
- `bargraph`/`linegraph`: Custom plotting functions (optional)

---

## Template Development

### Jinja2 Template Basics

**Main template file (templates/base.html):**

```html
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <title>{{ report.title }} - MultiQC Report</title>

    {# Include default CSS #} {% include 'header_css.html' %} {# Custom CSS #}
    <link rel="stylesheet" href="assets/css/custom.css" />
  </head>
  <body>
    {# Custom header #} {% include 'header.html' %} {# Main report content from default template #} {{
    report.plot_compressed_json | safe }}

    <div class="container-fluid main-content">
      {% for section in report.sections %}
      <div id="mqc-section-{{ section['id'] }}" class="mqc-section">
        <h2>{{ section['name'] }}</h2>
        <div class="mqc-section-content">{{ section['content'] | safe }}</div>
      </div>
      {% endfor %}
    </div>

    {# Custom footer #} {% include 'footer.html' %} {# Default JavaScript #} {% include 'footer_javascript.html' %} {#
    Custom JavaScript #}
    <script src="assets/js/custom.js"></script>
  </body>
</html>
```

### Available Template Variables

**report namespace:**

```python
report.title          # Report title
report.intro          # Introduction text
report.sections       # List of report sections
report.general_stats  # General statistics table
report.data_sources   # List of input files
report.plot_data      # Plot data structures
report.plot_compressed_json  # Compressed plot data for JavaScript
```

**config namespace:**

```python
config.title           # Configured report title
config.subtitle        # Report subtitle
config.intro_text      # Custom introduction
config.report_comment  # Additional comments
config.data_format     # Output data format
config.plots_force_flat  # Force static plots
```

**section structure:**

```python
section = {
    'id': 'fastqc',           # Section identifier
    'name': 'FastQC',         # Display name
    'content': '<div>...</div>',  # HTML content
    'plot_data': {...},       # Associated plot data
    'description': '...',     # Section description
}
```

### Child Template Customization

**Custom header (templates/header.html):**

```html
<header class="navbar navbar-default navbar-fixed-top">
  <div class="container-fluid">
    <div class="navbar-header">
      {# Institutional logo #}
      <a class="navbar-brand" href="https://your-institution.org">
        <img src="assets/img/logo.png" alt="Institution Logo" style="height: 40px; margin-top: -10px;" />
      </a>

      {# Report title #}
      <span class="navbar-text"> {{ report.title }} </span>
    </div>

    {# Additional navigation items #}
    <div class="navbar-right">
      <span class="navbar-text"> Generated: {{ report.creation_date }} </span>
    </div>
  </div>
</header>
```

**Custom footer (templates/footer.html):**

```html
<footer class="footer">
  <div class="container-fluid">
    <p class="text-muted text-center">
      Report generated by
      <a href="https://multiqc.info">MultiQC</a>
      using {{ config.template }} template.
      <br />
      <a href="https://your-institution.org">Your Institution</a> | Contact:
      <a href="mailto:support@your-institution.org">support@your-institution.org</a>
    </p>
  </div>
</footer>
```

### Helper Functions

**Include external files:**

```html
{# Include file contents directly #} {% include 'header.html' %} {# Include file with base64 encoding (for images in
CSS) #} {{ include_file('assets/img/background.png', b64=True) }}
```

**Material icons (default template):**

```html
{# Add SVG icon #} {{ material_icon('check_circle', size=24, color='#28a745') }} {{ material_icon('warning', size=20,
color='#ffc107') }} {{ material_icon('error', size=20, color='#dc3545') }}
```

**Common icon names:**

- `check_circle`: Success/complete
- `warning`: Warning/caution
- `error`: Error/failure
- `info`: Information
- `help`: Help/documentation
- `settings`: Configuration
- `download`: Download action

---

## Asset Management

### CSS Customization

**Custom styles (assets/css/custom.css):**

```css
/* Institutional color scheme */
:root {
  --primary-color: #0066cc;
  --secondary-color: #6c757d;
  --success-color: #28a745;
  --warning-color: #ffc107;
  --danger-color: #dc3545;
  --institution-font: "Open Sans", sans-serif;
}

/* Custom header styling */
.navbar-default {
  background-color: var(--primary-color);
  border-color: var(--primary-color);
}

.navbar-default .navbar-brand,
.navbar-default .navbar-text {
  color: white !important;
}

/* Custom section headers */
.mqc-section h2 {
  color: var(--primary-color);
  border-bottom: 2px solid var(--primary-color);
  padding-bottom: 10px;
  margin-bottom: 20px;
  font-family: var(--institution-font);
}

/* Sample name highlighting */
.sample-name {
  font-weight: bold;
  color: var(--primary-color);
}

/* Custom table styling */
.table-custom {
  border: 1px solid var(--secondary-color);
}

.table-custom thead {
  background-color: var(--primary-color);
  color: white;
}

/* Responsive adjustments */
@media (max-width: 768px) {
  .navbar-brand img {
    height: 30px !important;
  }
}
```

### JavaScript Customization

**Custom interactions (assets/js/custom.js):**

```javascript
// Wait for MultiQC to initialize
$(document).ready(function () {
  console.log("Custom template loaded");

  // Add custom click handlers
  $(".sample-name").on("click", function () {
    var sampleId = $(this).data("sample-id");
    console.log("Clicked sample:", sampleId);
    // Custom behavior here
  });

  // Customize plot options
  if (typeof mqc_plots !== "undefined") {
    // Modify plot configurations
    $.each(mqc_plots, function (plotId, plotConfig) {
      // Add institutional branding to plots
      if (plotConfig.config && plotConfig.config.layout) {
        plotConfig.config.layout.font = {
          family: "Open Sans, sans-serif",
        };
      }
    });
  }

  // Add custom toolbar buttons
  $("#mqc-toolbar").append(
    '<button class="btn btn-default btn-sm" onclick="exportCustomFormat()">' + "Export Custom Format</button>",
  );
});

// Custom export function
function exportCustomFormat() {
  // Implementation for custom data export
  alert("Exporting custom format...");
}
```

### Image Assets

**Logo preparation:**

- **Format**: PNG with transparency (for flexibility)
- **Size**: 2x resolution for retina displays (e.g., 400x80 for 200x40 display)
- **Location**: `assets/img/logo.png`
- **Usage**: `<img src="assets/img/logo.png" alt="Logo">`

**Background images:**

- Optimize file size (compress)
- Use WebP format for modern browsers
- Provide PNG fallback

**Embedding in CSS:**

```css
.custom-background {
  background-image: url("../img/background.png");
  background-size: cover;
}
```

---

## Custom Plotting

### Override Default Plot Functions

**Bar graph customization:**

```python
# In __init__.py or separate module

def bargraph(plotdata, pconfig):
    """
    Custom bar graph implementation using different library.

    Args:
        plotdata: Plot data structure from MultiQC
        pconfig: Plot configuration dictionary

    Returns:
        HTML string with custom plot, or None for default
    """
    # Example: Use alternative plotting library
    try:
        import altair as alt
        import pandas as pd

        # Convert plotdata to DataFrame
        df = pd.DataFrame(plotdata)

        # Create Altair chart
        chart = alt.Chart(df).mark_bar().encode(
            x='sample:N',
            y='value:Q',
            color='category:N',
            tooltip=['sample', 'value', 'category']
        ).properties(
            width=600,
            height=400,
            title=pconfig.get('title', 'Custom Bar Graph')
        )

        # Return HTML with embedded chart
        return chart.to_html()

    except ImportError:
        # Fallback to default if altair not available
        return None
```

**Line graph customization:**

```python
def linegraph(plotdata, pconfig):
    """
    Custom line graph with institutional styling.
    """
    # Return None to use default with your CSS styling
    return None
```

### Plot Configuration

**Modify plot appearance:**

```python
# In your custom template module

from multiqc import config

# Force all plots to be static images (no interactivity)
config.plots_force_flat = True

# Or set via MultiQC config file (multiqc_config.yaml):
# plots_force_flat: True
# plots_force_interactive: False
```

---

## Development Workflow

### Development Mode

**Enable rapid iteration:**

```bash
# Install template in development mode
cd multiqc_<template-name>
pip install -e .

# Run MultiQC with development flag
multiqc --development data/ --template <template-name>
```

**Development mode benefits:**

- Loads CSS/JS directly from source (no rebuild needed)
- Links plot images externally
- Exports uncompressed JSON plot data
- Immediate feedback on template changes

**Development cycle:**

1. **Edit template files** (base.html, CSS, JS)
2. **Run MultiQC** with `--development` flag
3. **Refresh browser** to see changes (no reinstall needed)
4. **Iterate** until satisfied
5. **Test production build** without `--development` flag

### Local Testing

**Test with sample data:**

```bash
# Clone MultiQC for test data
git clone https://github.com/ewels/MultiQC.git
cd MultiQC/test_data

# Run with your template
multiqc . --template <template-name> --development

# Output: multiqc_report.html
```

**Test different configurations:**

```bash
# Test with custom config
multiqc data/ --template <template-name> \
  --config multiqc_config.yaml \
  --development

# Test with different data
multiqc other_data/ --template <template-name> --development
```

### Validation Checklist

- [ ] Template loads without errors
- [ ] All sections render correctly
- [ ] Plots display properly (both interactive and static)
- [ ] Custom CSS applies correctly
- [ ] Custom JavaScript executes without errors
- [ ] Logo/images display at correct size
- [ ] Responsive design works (test mobile view)
- [ ] All links work correctly
- [ ] Footer information is accurate
- [ ] Export functions work (if implemented)
- [ ] Browser console shows no errors
- [ ] Template works with various data types
- [ ] Child template inherits parent functions correctly

---

## Testing and Validation

### Test Suite

**Create test script (test_template.sh):**

```bash
#!/bin/bash
# Test custom MultiQC template

TEMPLATE_NAME="<template-name>"
TEST_DATA="test_data"

echo "Testing MultiQC template: $TEMPLATE_NAME"

# Test 1: Basic rendering
echo "Test 1: Basic rendering..."
multiqc $TEST_DATA --template $TEMPLATE_NAME -f -o test_output/basic
if [ $? -eq 0 ]; then
    echo "✓ Basic rendering passed"
else
    echo "✗ Basic rendering failed"
    exit 1
fi

# Test 2: Development mode
echo "Test 2: Development mode..."
multiqc $TEST_DATA --template $TEMPLATE_NAME --development -f -o test_output/dev
if [ $? -eq 0 ]; then
    echo "✓ Development mode passed"
else
    echo "✗ Development mode failed"
    exit 1
fi

# Test 3: Custom config
echo "Test 3: Custom configuration..."
multiqc $TEST_DATA --template $TEMPLATE_NAME \
  --config test_config.yaml -f -o test_output/config
if [ $? -eq 0 ]; then
    echo "✓ Custom config passed"
else
    echo "✗ Custom config failed"
    exit 1
fi

# Test 4: Various data types
for data_dir in $TEST_DATA/fastqc $TEST_DATA/star $TEST_DATA/salmon; do
    if [ -d "$data_dir" ]; then
        echo "Test 4: Data type $(basename $data_dir)..."
        multiqc $data_dir --template $TEMPLATE_NAME -f \
          -o test_output/$(basename $data_dir)
        if [ $? -eq 0 ]; then
            echo "✓ $(basename $data_dir) passed"
        else
            echo "✗ $(basename $data_dir) failed"
            exit 1
        fi
    fi
done

echo ""
echo "All tests passed! ✓"
echo "Check test_output/ for generated reports"
```

### Visual Verification

**Checklist for manual review:**

1. **Open generated report in browsers:**

   - [ ] Chrome/Chromium
   - [ ] Firefox
   - [ ] Safari
   - [ ] Edge

2. **Check responsive design:**

   - [ ] Desktop (1920x1080)
   - [ ] Tablet (768x1024)
   - [ ] Mobile (375x667)

3. **Verify visual elements:**

   - [ ] Logo displays correctly
   - [ ] Colors match branding
   - [ ] Fonts load properly
   - [ ] Icons render correctly
   - [ ] Spacing is consistent

4. **Test interactivity:**

   - [ ] Plot tooltips work
   - [ ] Plot zooming/panning works
   - [ ] Sample filtering works
   - [ ] Custom buttons function
   - [ ] Links navigate correctly

5. **Check browser console:**
   - [ ] No JavaScript errors
   - [ ] No CSS warnings
   - [ ] No missing resources (404s)

### Common Issues and Solutions

**Issue: Template not found**

```bash
# Check installation
pip show multiqc_<template-name>

# Verify entry point
multiqc --list-templates
# Should show your template name

# Reinstall if needed
pip install -e . --force-reinstall
```

**Issue: Assets not loading**

```python
# Verify copy_files in __init__.py
copy_files = [
    'assets',  # Ensure this is present
]

# Check include_package_data in setup.py
include_package_data=True
```

**Issue: Child template not inheriting**

```python
# Verify template_parent in __init__.py
template_parent = 'default'  # Must match existing template

# Check parent template exists
multiqc --list-templates
```

**Issue: JavaScript errors**

```javascript
// Wrap in document ready
$(document).ready(function () {
  // Your code here - ensures MultiQC is initialized
});
```

---

## Distribution

### Package for Distribution

**Build package:**

```bash
# Update version in setup.py
# Commit all changes

# Build distribution files
python setup.py sdist bdist_wheel

# Check dist/ directory
ls dist/
# multiqc_<template-name>-1.0.0-py3-none-any.whl
# multiqc_<template-name>-1.0.0.tar.gz
```

### Installation Methods

**Method 1: PyPI (public distribution)**

```bash
# Upload to PyPI
pip install twine
twine upload dist/*

# Users install with:
pip install multiqc_<template-name>
multiqc data/ --template <template-name>
```

**Method 2: GitHub (open source)**

```bash
# Users install directly from GitHub:
pip install git+https://github.com/username/multiqc_<template-name>.git

# Or clone and install:
git clone https://github.com/username/multiqc_<template-name>.git
cd multiqc_<template-name>
pip install .
```

**Method 3: Private distribution**

```bash
# Share .whl file internally
# Users install with:
pip install multiqc_<template-name>-1.0.0-py3-none-any.whl
```

**Method 4: Direct template directory**

```bash
# Users can place template in:
# ~/.config/multiqc/templates/<template-name>/

# Or use --template-dir flag:
multiqc data/ --template <template-name> \
  --template-dir /path/to/templates/
```

### Documentation

**README.md template:**

````markdown
# MultiQC Template: <template-name>

Custom MultiQC template for [organization/purpose].

## Features

- Custom institutional branding
- Modified color scheme
- Enhanced visualizations
- [Other features]

## Installation

```bash
pip install multiqc_<template-name>
```
````

## Usage

```bash
multiqc data/ --template <template-name>
```

## Configuration

Optional configuration in `multiqc_config.yaml`:

```yaml
template: <template-name>

# Custom settings
custom_logo: true
color_scheme: blue
```

## Screenshots

[Include screenshots of the custom template]

## Development

```bash
git clone https://github.com/username/multiqc_<template-name>.git
cd multiqc_<template-name>
pip install -e .
multiqc test_data/ --template <template-name> --development
```

## License

[Your license]

## Credits

Built on [MultiQC](https://multiqc.info) by Phil Ewels.

````

### Contribution to MultiQC

**For general-purpose templates:**

1. Fork MultiQC repository
2. Add template to `multiqc/templates/<template-name>/`
3. Follow MultiQC contribution guidelines
4. Submit pull request
5. Benefits: Bundled with MultiQC, wider adoption

**For specialized templates:**

- Distribute as separate package
- Better for institutional/specific use cases
- Easier to maintain independently
- Can have different release cycle

---

## Complete Example: Lab Branding Template

### Project Setup

```bash
# Create project structure
mkdir -p multiqc_mylab/{multiqc_mylab/templates/assets/{css,js,img},test_data}
cd multiqc_mylab
````

### setup.py

```python
from setuptools import setup, find_packages

setup(
    name='multiqc_mylab',
    version='1.0.0',
    author='Lab Director',
    author_email='director@mylab.edu',
    description='Custom MultiQC template with MyLab branding',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/mylab/multiqc_mylab',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['multiqc>=1.14'],
    entry_points={
        'multiqc.templates.v1': [
            'mylab = multiqc_mylab'
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
)
```

### multiqc_mylab/**init**.py

```python
import os

template_parent = 'default'
template_dir = os.path.dirname(__file__)
base_fn = 'base.html'
copy_files = ['assets']
```

### multiqc_mylab/templates/base.html

```html
<!DOCTYPE html>
<html>
  <head>
    <title>{{ report.title }} - MyLab Report</title>
    {% include 'header_css.html' %}
    <link rel="stylesheet" href="assets/css/custom.css" />
  </head>
  <body>
    {% include 'header.html' %} {{ report.plot_compressed_json | safe }}

    <div class="container-fluid main-content">
      {% for section in report.sections %}
      <div id="mqc-section-{{ section['id'] }}" class="mqc-section">
        <h2>{{ section['name'] }}</h2>
        {{ section['content'] | safe }}
      </div>
      {% endfor %}
    </div>

    {% include 'footer.html' %} {% include 'footer_javascript.html' %}
  </body>
</html>
```

### multiqc_mylab/templates/header.html

```html
<header class="navbar navbar-default navbar-fixed-top">
  <div class="container-fluid">
    <div class="navbar-header">
      <a href="https://mylab.edu">
        <img src="assets/img/logo.png" alt="MyLab" style="height: 40px;" />
      </a>
      <span class="navbar-text">{{ report.title }}</span>
    </div>
  </div>
</header>
```

### multiqc_mylab/templates/footer.html

```html
<footer class="footer">
  <div class="container-fluid text-center">
    <p>Generated by <a href="https://multiqc.info">MultiQC</a></p>
    <p><a href="https://mylab.edu">MyLab</a> | <a href="mailto:support@mylab.edu">Contact</a></p>
  </div>
</footer>
```

### multiqc_mylab/templates/assets/css/custom.css

```css
:root {
  --mylab-blue: #003087;
  --mylab-gold: #ffc72c;
}

.navbar-default {
  background-color: var(--mylab-blue);
}

.navbar-text {
  color: white !important;
  margin-left: 15px;
}

.mqc-section h2 {
  color: var(--mylab-blue);
  border-bottom: 3px solid var(--mylab-gold);
}
```

### Development and Testing

```bash
# Install in development mode
pip install -e .

# Test with sample data
multiqc test_data/ --template mylab --development

# Production build test
multiqc test_data/ --template mylab -f

# Build distribution
python setup.py sdist bdist_wheel
```

---

## Best Practices

### Template Development

1. **Start with child template**: Inherit from `default` unless you need full control
2. **Use version control**: Git repository from the start
3. **Test iteratively**: Use `--development` flag for rapid iteration
4. **Document everything**: Clear README with examples
5. **Minimize assets**: Optimize images and CSS for faster loading
6. **Validate HTML/CSS**: Use validators to ensure standard compliance
7. **Test responsive design**: Ensure reports work on all screen sizes
8. **Follow MultiQC conventions**: Stay consistent with default template patterns

### Performance Optimization

1. **Compress images**: Use tools like ImageOptim, TinyPNG
2. **Minimize CSS/JS**: Use minification for production
3. **Lazy load images**: Defer non-critical images
4. **Use CSS sprites**: Combine small images
5. **Enable caching**: Set appropriate cache headers
6. **Avoid external dependencies**: Bundle resources when possible

### Maintenance

1. **Test with MultiQC updates**: Ensure compatibility with new releases
2. **Monitor GitHub issues**: Address user feedback
3. **Version appropriately**: Use semantic versioning
4. **Maintain changelog**: Document changes between versions
5. **Update dependencies**: Keep up with MultiQC version requirements

---

## Resources

### Official Documentation

- **MultiQC Docs**: https://docs.seqera.io/multiqc/
- **Template Development**: https://docs.seqera.io/multiqc/development/templates
- **GitHub Repository**: https://github.com/ewels/MultiQC

### Tutorials and Guides

- **Bytesize Talk**: https://nf-co.re/events/2022/bytesize-36-multiqc
- **Example Templates**: Browse `multiqc/templates/` in MultiQC repository

### Related Skills

- For finding MultiQC source: Use `investigating-bioinformatics-tool` skill
- For analyzing MultiQC code: Use `analyzing-tool-implementation` skill
- For MultiQC documentation: Use `searching-tool-online-docs` skill

---

## Quick Reference

### Template Creation Workflow

```bash
# 1. Create project structure
mkdir -p multiqc_<name>/{multiqc_<name>/templates/assets/{css,js,img},test_data}

# 2. Create setup.py with entry point
# 3. Create __init__.py with template configuration
# 4. Create templates/base.html (or custom components)
# 5. Add custom CSS/JS/images in assets/

# 6. Install in development mode
cd multiqc_<name>
pip install -e .

# 7. Test iteratively
multiqc test_data/ --template <name> --development

# 8. Validate and refine
# 9. Build distribution
python setup.py sdist bdist_wheel

# 10. Distribute (PyPI, GitHub, or internal)
```

### Common Commands

```bash
# List available templates
multiqc --list-templates

# Use custom template
multiqc data/ --template <template-name>

# Development mode
multiqc data/ --template <template-name> --development

# Custom template directory
multiqc data/ --template-dir /path/to/templates/ --template <name>

# With custom config
multiqc data/ --template <name> --config multiqc_config.yaml
```

---

This skill provides comprehensive guidance for creating custom MultiQC templates. Use the development workflow with validation checkpoints to ensure your template works correctly before distribution.
