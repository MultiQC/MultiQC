FROM python:3.11-slim

LABEL author="Phil Ewels" \
      description="MultiQC" \
      maintainer="phil.ewels@seqera.io"

# Add the MultiQC source files to the container
ADD . /usr/src/multiqc
WORKDIR /usr/src/multiqc

# - Install `ps` for Nextflow
# - Install MultiQC through pip
# - Delete unnecessary Python files
# - Remove MultiQC source directory
# - Add custom group and user
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/* && \
    pip install --upgrade pip && \
    pip install -v --no-cache-dir . && \
    apt-get clean -y && \
    for cache_dir in $(find /usr/local/lib/python3.11 \( -iname '*.c' -o -iname '*.pxd' -o -iname '*.pyd' -o -iname '__pycache__' \)); do \
      rm -rf "$cache_dir"; \
    done; \
    rm -rf "/usr/src/multiqc/" && \
    groupadd --gid 1000 multiqc && \
    useradd -ms /bin/bash --create-home --gid multiqc --uid 1000 multiqc

# Set to be the new user
USER multiqc

# Set default workdir to user home
WORKDIR /home/multiqc