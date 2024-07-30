#!/usr/bin/env bash

echo "Previewing docs site locally..."
docker run -i -t -p 3000:3000 -v ${PWD}:/MultiQC ghcr.io/multiqc/website:latest
