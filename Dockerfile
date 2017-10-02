FROM python:2.7-slim

LABEL \
  author="Phil Ewels" \
  description="MultiQC" \
  maintainer="phil.ewels@scilifelab.se"

# Setup ENV variables
ENV MULTIQC_VERSION="1.1"

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
  g++ \
  git \
  wget \
  && wget --quiet -O /opt/get-pip.py https://bootstrap.pypa.io/get-pip.py \
  && python /opt/get-pip.py \
  && rm -rf /var/lib/apt/lists/* /opt/get-pip.py

# Install MultiQC
RUN \
  pip install networkx==1.11 \
  && pip install git+git://github.com/ewels/MultiQC.git@v${MULTIQC_VERSION}
