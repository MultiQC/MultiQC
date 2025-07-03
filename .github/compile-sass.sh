#!/bin/bash
# Pre-commit hook to compile Sass when SCSS files are modified

# Check if any SCSS files in the default template were modified
if git diff --cached --name-only | grep -E "^multiqc/templates/default/scss/.*\.scss$"; then
    echo "SCSS files modified, compiling CSS..."

    # Change to the default template directory
    cd multiqc/templates/default

    # Check if node_modules exists
    if [ ! -d "node_modules" ]; then
        echo "Installing npm dependencies..."
        npm install
    fi

    # Build the CSS
    npm run build-css

    # Add the compiled CSS to the commit
    git add assets/css/multiqc.min.css

    echo "CSS compilation complete."
fi
