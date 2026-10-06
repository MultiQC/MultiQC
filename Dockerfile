FROM python:3.13-slim

LABEL author="Phil Ewels & Vlad Savelyev" \
      description="MultiQC" \
      maintainer="phil.ewels@seqera.io"

# Optional pandoc installation for PDF support
ARG INSTALL_PANDOC=false

# - Install `ps` for Nextflow
# - Install MultiQC through pip
# - Delete unnecessary Python files
# - Remove build artifacts
# - Add custom group and user
# Source is bind-mounted rather than COPY'd so it never ends up in an image layer
RUN \
    --mount=type=bind,source=multiqc,target=/usr/src/multiqc/multiqc \
    --mount=type=bind,source=pyproject.toml,target=/usr/src/multiqc/pyproject.toml \
    --mount=type=bind,source=README.md,target=/usr/src/multiqc/README.md \
    --mount=type=bind,source=LICENSE,target=/usr/src/multiqc/LICENSE \
    echo "Docker build log: Run apt-get update" 1>&2 && \
    apt-get update -y -qq \
    && \
    echo "Docker build log: Install procps" 1>&2 && \
    apt-get install -y -qq procps && \
    if [ "$INSTALL_PANDOC" = "true" ]; then \
        echo "Docker build log: Install pandoc and LaTeX for PDF generation" 1>&2 && \
        apt-get install -y -qq pandoc texlive-latex-base texlive-fonts-recommended texlive-latex-extra texlive-luatex; \
    fi && \
    echo "Docker build log: Clean apt cache" 1>&2 && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean -y && \
    echo "Docker build log: Upgrade pip and install multiqc" 1>&2 && \
    pip install --quiet --upgrade pip && \
    #################
    # Install MultiQC
    pip install --verbose --no-cache-dir /usr/src/multiqc && \
    echo "Docker build log: Delete python cache directories" 1>&2 && \
    find /usr/local/lib/python3.13 \( -iname '*.c' -o -iname '*.pxd' -o -iname '*.pyd' -o -iname '__pycache__' \) -printf "\"%p\" " | \
    xargs rm -rf {} && \
    echo "Docker build log: Delete build artifacts" 1>&2 && \
    rm -rf /usr/src/multiqc/build /usr/src/multiqc/*.egg-info && \
    echo "Docker build log: Add multiqc user and group" 1>&2 && \
    groupadd --gid 1000 multiqc && \
    useradd -ms /bin/bash --create-home --gid multiqc --uid 1000 multiqc

# Set to be the new user
USER multiqc

# Set default workdir to user home
WORKDIR /home/multiqc

# Check everything is working smoothly
RUN echo "Docker build log: Testing multiqc" 1>&2 && \
    multiqc --help

# Display the command line help if the container is run without any parameters
CMD multiqc --help
