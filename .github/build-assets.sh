#!/bin/bash
# Pre-commit hook to build assets with Vite when source files are modified

# Check if any source files in the default template were modified
if git diff --cached --name-only | grep -E "^multiqc/templates/default/src/.*\.(scss|js)$"; then
    echo "Source files modified, building assets with Vite..."

    # Change to the default template directory
    cd multiqc/templates/default

    # Check if node_modules exists
    if [ ! -d "node_modules" ]; then
        echo "Installing npm dependencies..."
        npm install
    fi

    # Build with Vite
    npm run build

    echo "Build complete."
fi
