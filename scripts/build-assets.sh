#!/bin/bash
# Pre-commit hook to build assets with Vite when source files are modified

# Check if any source files were passed as arguments (pre-commit filters by the files pattern)
if [ $# -gt 0 ]; then
    echo "Source files modified: $*"
    echo "Building assets with Vite..."

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
else
    echo "No relevant source files modified."
fi
