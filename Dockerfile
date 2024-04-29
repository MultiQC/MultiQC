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
RUN \
    echo "Docker build log: Run apt-get update" 1>&2 && \
    apt-get update -y -qq \
    && \
    echo "Docker build log: Install procps" 1>&2 && \
    apt-get install -y -qq procps && \
    echo "Docker build log: Clean apt cache" 1>&2 && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean -y && \
    echo "Docker build log: Upgrade pip and install multiqc" 1>&2 && \
    pip install --quiet --upgrade pip && \
    pip install --verbose --no-cache-dir . && \
    echo "Docker build log: Delete python cache directories" 1>&2 && \
    find /usr/local/lib/python3.11 \( -iname '*.c' -o -iname '*.pxd' -o -iname '*.pyd' -o -iname '__pycache__' \) -printf "\"%p\" " | \
    xargs rm -rf {} && \
    echo "Docker build log: Delete /usr/src/multiqc" 1>&2 && \
    rm -rf "/usr/src/multiqc/" && \
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
