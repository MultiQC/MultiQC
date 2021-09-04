FROM alpine:3.13

# Based on teeny-tiny matplotlib image based on alpine (czentye/matplotlib-minimal:3.1.2)
# But rebuilding here since alpine:latest is a multi-arch build

LABEL author="Phil Ewels" \
      description="MultiQC" \
      maintainer="phil.ewels@scilifelab.se"

ARG MATPLOTLIB_VERSION="3.1.2"

ARG USER_NAME="multiqc_user"
ARG USER_ID="1000"
ARG GROUP_NAME="multiqc_user"
ARG GROUP_ID="1000"

# Add the MultiQC source files to the container
ADD . /usr/src/multiqc
WORKDIR /usr/src/multiqc

RUN apk add --repository http://dl-cdn.alpinelinux.org/alpine/edge/main \
            --update \
            --no-cache \
            bash \
            build-base \
            freetype \
            freetype-dev \
            libgfortran \
            libpng \
            libpng-dev \
            libstdc++ \
            python3 \
            python3-dev \
            py3-setuptools \
            py3-pip && \
    ln -fs /usr/include/locale.h /usr/include/xlocale.h && \
    ln -fs /usr/bin/python3 /usr/local/bin/python && \
    ln -fs /usr/bin/pip3 /usr/local/bin/pip && \
    pip3 install -v --no-cache-dir --upgrade \
      pip && \
    pip3 install -v --no-cache-dir \
      numpy && \
    pip3 install -v --no-cache-dir \
      matplotlib=="$MATPLOTLIB_VERSION" && \
    pip3 install -v --no-cache-dir \
      wheel && \
    pip3 install -v --no-cache-dir \
      . && \
    apk del \
      --purge \
      build-base \
      libgfortran \
      libpng-dev \
      freetype-dev \
      python3-dev && \
    rm -vrf /var/cache/apk/* && \
    addgroup \
      -g "${GROUP_ID}" \
      "${GROUP_NAME}" && \
    adduser \
      -D \
      -G "${GROUP_NAME}" \
      -u "${USER_ID}" \
      "${USER_NAME}"

# Add the MultiQC source files to the container
USER "${USER_NAME}"
  
# Set up entrypoint and cmd for easy docker usage
CMD [ "multiqc" ]

