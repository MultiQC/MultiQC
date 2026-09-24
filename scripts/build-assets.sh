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

    # The 'disco' template re-exports the default template's JS, so its own bundle
    # goes stale whenever the default sources change.
    if [ -d "../disco" ]; then
        cd ../disco
        if [ ! -d "node_modules" ]; then
            npm ci
        fi
        npm run build
        cd ../default
    fi

    # 'original' loads plain script tags rather than a bundle, so it needs its own copy
    # of DOMPurify. Copy it from node_modules so package.json stays the only version pin.
    cp node_modules/dompurify/dist/purify.min.js ../original/assets/js/packages/dompurify.min.js

    echo "Build complete."
else
    echo "No relevant source files modified."
fi
