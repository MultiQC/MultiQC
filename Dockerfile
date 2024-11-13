# Build stage
FROM python:3.12-slim as builder

# Add the MultiQC source files to the container
RUN mkdir /usr/src/multiqc
COPY LICENSE /usr/src/multiqc/
COPY README.md /usr/src/multiqc/
COPY multiqc /usr/src/multiqc/multiqc
COPY pyproject.toml /usr/src/multiqc/
COPY MANIFEST.in /usr/src/multiqc/
COPY scripts /usr/src/multiqc/scripts
COPY setup.py /usr/src/multiqc/
COPY tests /usr/src/multiqc/tests

WORKDIR /usr/src/multiqc
RUN pip install --no-cache-dir .

# Runtime stage
FROM python:3.12-slim
LABEL author="Phil Ewels & Vlad Savelyev" \
      description="MultiQC" \
      maintainer="phil.ewels@seqera.io"

COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin/multiqc /usr/local/bin/multiqc

RUN groupadd --gid 1000 multiqc && \
    useradd -ms /bin/bash --create-home --gid multiqc --uid 1000 multiqc

# Set to be the new user
USER multiqc

# Set default workdir to user home
WORKDIR /home/multiqc

# Check everything is working smoothly
RUN echo "Docker build log: Testing multiqc" 1>&2 && \
    multiqc --help 

# Display the command line help if the container is run without any parameters
CMD ["multiqc", "--help"]
