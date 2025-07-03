# MultiQC Default Template

This directory contains the default HTML template for MultiQC reports. The template uses Bootstrap 5.3.7 and custom styles compiled from Sass.

## Directory Structure

```
default/
├── scss/
│   ├── main.scss      # Main file that imports Bootstrap and custom styles
│   └── custom.scss    # Custom MultiQC styles
├── assets/
│   ├── css/
│   │   └── multiqc.min.css  # Compiled CSS (do not edit directly)
│   ├── js/
│   └── img/
├── node_modules/      # npm dependencies (not in git)
└── package.json       # npm configuration
```

## Development Setup

1. Install dependencies:

   ```bash
   cd multiqc/templates/default
   npm install
   ```

2. Build CSS from Sass:

   ```bash
   npm run build-css
   ```

3. Watch for changes during development:
   ```bash
   npm run watch-css
   ```

## Making Style Changes

1. Edit the Sass files in the `scss/` directory:

   - `custom.scss` - Contains all custom MultiQC styles
   - `main.scss` - Imports Bootstrap and custom styles

2. The CSS will be automatically compiled when:
   - You run `npm run build-css` manually
   - You have `npm run watch-css` running
   - You commit changes (via pre-commit hook)

## Important Notes

- **Never edit `multiqc.min.css` directly** - it is automatically generated
- The compiled CSS is minified and includes both Bootstrap 5.3.7 and custom styles
- The pre-commit hook ensures the compiled CSS is always up-to-date in the repository
- All assets are base64-encoded into the HTML report to make it standalone
