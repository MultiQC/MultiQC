ARG PYTHON_IMAGE=python:3.13-slim

FROM ${PYTHON_IMAGE} AS builder

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /usr/src/multiqc

# Only copy sources needed to build the package wheel.
COPY LICENSE README.md pyproject.toml MANIFEST.in setup.py ./
COPY docs ./docs
COPY multiqc ./multiqc
COPY scripts ./scripts

RUN python -m pip install --upgrade pip setuptools wheel \
    && python -m pip wheel --no-deps --wheel-dir /tmp/wheels .


FROM ${PYTHON_IMAGE}

LABEL org.opencontainers.image.title="MultiQC" \
      org.opencontainers.image.description="Summarize analysis results from multiple tools and samples in a single report" \
      org.opencontainers.image.authors="Phil Ewels, Vlad Savelyev" \
      org.opencontainers.image.source="https://github.com/MultiQC/MultiQC"

# Optional pandoc installation for PDF support
ARG INSTALL_PANDOC=false

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# - Install runtime dependencies only (procps for Nextflow, optional pandoc+latex for PDF)
# - Use no-install-recommends to minimize image size and attack surface
# - Clean apt lists in the same layer
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && if [ "$INSTALL_PANDOC" = "true" ]; then \
        apt-get install -y --no-install-recommends \
            pandoc \
            texlive-latex-base \
            texlive-fonts-recommended \
            texlive-latex-extra \
            texlive-luatex; \
    fi \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /tmp/wheels /tmp/wheels

RUN python -m pip install --upgrade pip \
    && python -m pip install /tmp/wheels/*.whl \
    && rm -rf /tmp/wheels

# Run as unprivileged user by default.
RUN groupadd --gid 1000 multiqc \
    && useradd --uid 1000 --gid multiqc --create-home --home-dir /home/multiqc --shell /usr/sbin/nologin multiqc

USER 1000:1000

WORKDIR /home/multiqc

# Check everything is working smoothly.
RUN multiqc --help > /dev/null

# Display command line help if the container is run without parameters.
CMD ["multiqc", "--help"]
