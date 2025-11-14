FROM python:3.13-slim

LABEL author="Phil Ewels & Vlad Savelyev" \
      description="MultiQC" \
      maintainer="phil.ewels@seqera.io"

# Optional pandoc installation for PDF support
ARG INSTALL_PANDOC=false

# Copy uv from the official image for faster dependency resolution
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Install system dependencies
RUN echo "Docker build log: Run apt-get update" 1>&2 && \
    apt-get update -y -qq && \
    echo "Docker build log: Install procps" 1>&2 && \
    apt-get install -y -qq procps && \
    if [ "$INSTALL_PANDOC" = "true" ]; then \
        echo "Docker build log: Install pandoc and LaTeX for PDF generation" 1>&2 && \
        apt-get install -y -qq pandoc texlive-latex-base texlive-fonts-recommended texlive-latex-extra texlive-luatex; \
    fi && \
    echo "Docker build log: Clean apt cache" 1>&2 && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean -y

# Set up working directory
RUN mkdir /usr/src/multiqc
WORKDIR /usr/src/multiqc

# Copy dependency files first for better layer caching
COPY pyproject.toml uv.lock ./

# Install dependencies using uv (cached layer if dependencies don't change)
RUN echo "Docker build log: Install dependencies with uv" 1>&2 && \
    uv pip install --system --no-cache -r pyproject.toml

# Copy the MultiQC source files to the container
COPY LICENSE README.md MANIFEST.in setup.py ./
COPY docs ./docs
COPY multiqc ./multiqc
COPY scripts ./scripts
COPY tests ./tests

# Install MultiQC package using uv
RUN echo "Docker build log: Install MultiQC with uv" 1>&2 && \
    uv pip install --system --no-cache .

# Clean up build artifacts and Python cache
RUN echo "Docker build log: Delete python cache directories" 1>&2 && \
    find /usr/local/lib/python3.13 \( -iname '*.c' -o -iname '*.pxd' -o -iname '*.pyd' -o -iname '__pycache__' \) -printf "\"%p\" " | \
    xargs rm -rf {} && \
    echo "Docker build log: Delete /usr/src/multiqc" 1>&2 && \
    rm -rf "/usr/src/multiqc/"

# Add custom group and user
RUN echo "Docker build log: Add multiqc user and group" 1>&2 && \
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
