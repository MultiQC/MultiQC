#!/usr/bin/env python
"""
MultiQC is a tool to aggregate bioinformatics results across many samples into
a single report. It is written in Python and contains modules for a number of
common tools (FastQC, Bowtie, Picard and many others).

You can install MultiQC from PyPI as follows:

    pip install multiqc

Then it's just a case of going to your analysis directory and running the script:

    multiqc .

MultiQC will scan the specified directory ('.' is the current dir)
and produce a report detailing whatever it finds.

The report is created in `multiqc_report/multiqc_report.html` by default.
A zip file of the report is also generated to facilitate easy transfer and sharing.

Tab-delimited data files are created in `multiqc_report/report_data/` to give
easy access for downstream processing.

For more detailed instructions, run `multiqc -h`
"""

from setuptools import setup, find_packages

version = '0.2.1dev'

setup(
    name = 'multiqc',
    version = version,
    author = 'Phil Ewels',
    author_email = 'phil.ewels@scilifelab.se',
    description = "Create aggregate bioinformatics analysis report across many samples",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/ewels/MultiQC',
    download_url = 'https://github.com/ewels/MultiQC/releases',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    zip_safe = False,
    scripts = ['scripts/multiqc'],
    install_requires = [
        'jinja2',
        'simplejson',
        'pyyaml'
    ],
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)

