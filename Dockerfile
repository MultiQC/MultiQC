# Teeny-tiny matplotlib image based on alpine
FROM czentye/matplotlib-minimal:3.0.3

LABEL author="Phil Ewels" \
      description="MultiQC" \
      maintainer="phil.ewels@scilifelab.se"


# Add the MultiQC source files to the container
ADD . /usr/src/multiqc
WORKDIR /usr/src/multiqc

# Remove matplotlib requirement needed for py2 support
# TODO: We can get rid of this when MultiQC is py3 only
RUN sed -i 's/matplotlib>=2.1.1,<3.0.0/matplotlib>=2.1.1/g' setup.py

# Install MultiQC
RUN python setup.py install

# Set up entrypoint and cmd for easy docker usage
ENTRYPOINT [ "multiqc" ]
CMD [ "." ]
