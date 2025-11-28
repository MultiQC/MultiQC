# CLAUDE.md - MultiQC Default Template

This file provides guidance for working with the MultiQC default template frontend assets.

## Overview

The default template is the main HTML report template for MultiQC. It uses Bootstrap 5.3.7 for styling and layout, with custom JavaScript for interactive features like plot rendering, sample filtering, and toolbox controls.

## Directory Structure

```
multiqc/templates/default/
├── assets/           # Third-party libraries (jQuery, Plotly.js, Bootstrap icons)
├── compiled/         # ⚠️ AUTO-GENERATED - Do not edit directly!
│   ├── css/         # Bundled CSS output from Vite
│   └── js/          # Bundled JS output from Vite
├── src/
│   ├── js/          # Source JavaScript files
│   └── scss/        # Source SCSS files
├── *.html           # Jinja2 template files
├── package.json     # Node dependencies and build scripts
├── vite.config.js   # Vite bundler configuration
└── tsconfig.json    # TypeScript configuration
```

## Build System

The template uses **Vite** to bundle JavaScript and SCSS files into the `compiled/` directory.

### Setup

```bash
cd multiqc/templates/default
npm install
```

### Build Commands

```bash
# One-time build (for production or testing)
npm run build

# Watch mode (auto-rebuild on file changes during development)
npm run watch
```

### Important Rules

- **NEVER edit files in `compiled/`** - they are auto-generated and will be overwritten
- **Always edit source files** in `src/js/` and `src/scss/`
- **NEVER add inline styles** like `style="--bs-btn-padding-y: .25rem;`, use `custom.scss` if you must.
- Run `npm run build` after making changes before committing
- The build is also run automatically by pre-commit hooks

## Source Files

### JavaScript (`src/js/`)

Main files:

- **`multiqc.js`**: Entry point, imports all other modules
- **`toolbox.js`**: Toolbar controls (export buttons, highlighting, filtering)
- **`plots.js`**: Plotly.js plot decompression and rendering
- **`table.js`**: DataTables initialization and configuration
- **`highlight.js`**: Sample name highlighting/filtering across plots and tables

### SCSS (`src/scss/`)

- **`multiqc.scss`**: Main entry point
- **`variables.scss`**: Bootstrap variable overrides (colors, spacing, etc.)
- **`custom.scss`**: Custom styles for MultiQC-specific components

The SCSS imports Bootstrap 5.3.7 from npm and applies custom variables and styles.

## Assets

Third-party libraries in `assets/` are kept separate and base64-encoded directly into the final HTML report:

- jQuery 3.7.1
- Plotly.js (basic bundle)
- Bootstrap Icons
- DataTables

These are NOT bundled by Vite to maintain version control and reduce build complexity.

## Testing Changes

After making changes to JavaScript or SCSS:

1. Run `npm run build` to generate new compiled assets
2. Run MultiQC on test data: `multiqc test_data/`
3. Open the generated `multiqc_report.html` in a browser
4. Test your changes across different browsers if possible

## Common Tasks

### Changing Bootstrap Variables

Edit `src/scss/variables.scss`:

```scss
// Example: Change primary color
$primary: #3498db;
```

Then rebuild: `npm run build`

### Adding New JavaScript Functionality

1. Edit or create file in `src/js/`
2. Import in `src/js/multiqc.js` if it's a new file
3. Rebuild: `npm run build`
4. Test the generated report

### Modifying Template HTML

Edit the `.html` Jinja2 template files directly:

- `head.html`: `<head>` content, script/style tags
- `nav.html`: Navigation sidebar
- `base.html`: Main layout wrapper
- `footer.html`: Report footer

No build step needed for HTML changes (they use Jinja2 templating).

## Development Workflow

Best workflow for active development:

```bash
# Terminal 1: Watch for changes to JS/CSS
cd multiqc/templates/default
npm run watch

# Terminal 2: Run MultiQC with --development flag
multiqc /path/to/test/data --development --force
```

The watch mode will automatically rebuild assets when you save changes to source files.

### Auto-regenerate on HTML Changes

For HTML template changes (which aren't watched by `npm run watch`), use watchmedo:

```bash
# Install watchdog
pip install watchdog

# Watch for HTML changes and regenerate report
watchmedo shell-command \
  --patterns="*.html" \
  --recursive \
  --command='multiqc /path/to/test/data --development --force' \
  multiqc/templates/default/
```

## Browser Compatibility

The template targets modern browsers (ES6+). Key features used:

- Arrow functions
- Template literals
- `const`/`let`
- `fetch` API

Legacy browser support (IE11) was dropped in MultiQC v1.15.
