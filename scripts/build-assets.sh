#!/bin/bash
# Pre-commit hook to build assets with Vite when source files are modified

# Check if any source files were passed as arguments (prek filters by the files pattern)
if [ $# -gt 0 ]; then
    echo "Source files modified: $*"
    echo "Building assets with Vite..."

    # Change to the default template directory
    cd multiqc/templates/default

    # Install from the lockfile if needed. Use `npm ci`, not `npm install`:
    # `npm install` rewrites package-lock.json (differently per platform / npm
    # version), which dirties the tree and fails the hook in CI. `npm ci`
    # installs exactly what the lockfile pins and never modifies it.
    if [ ! -d "node_modules" ]; then
        echo "Installing npm dependencies..."
        npm ci
    fi

    # Build with Vite
    npm run build

    echo "Build complete."
else
    echo "No relevant source files modified."
fi
