#!/bin/bash
# Pre-commit hook to build assets with Vite when source files are modified

# Check if any source files were passed as arguments (prek filters by the files pattern)
# Install from the lockfile when node_modules is missing or older than it. Use
# `npm ci`, not `npm install`: `npm install` rewrites package-lock.json
# (differently per platform / npm version), which dirties the tree and fails the
# hook in CI. A stale install would build different output from CI.
install_deps() {
    if [ ! -f node_modules/.package-lock.json ] || [ package-lock.json -nt node_modules/.package-lock.json ]; then
        echo "Installing npm dependencies..."
        npm ci
    fi
}

if [ $# -gt 0 ]; then
    echo "Source files modified: $*"
    echo "Building assets with Vite..."

    # Change to the default template directory
    cd multiqc/templates/default

    install_deps

    # Build with Vite
    npm run build

    # The 'disco' template re-exports the default template's JS, so its own bundle
    # goes stale whenever the default sources change.
    if [ -d "../disco" ]; then
        cd ../disco
        install_deps
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
