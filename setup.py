#!/usr/bin/env python
"""
MultiQC is a tool to aggregate bioinformatics results across many samples into a single report. It is written in Python and contains modules for a number of common tools (FastQC, Bowtie, Picard and many others).

You can install MultiQC from PyPI as follows::

    pip install multiqc

Then it's just a case of going to your analysis directory and running the script::

    multiqc .

MultiQC will scan the specified directory (:code:`'.'` is the current dir) and produce a report detailing whatever it finds.

The report is created in :code:`multiqc_report.html` by default. Tab-delimited data files are created in :code:`multiqc_data/` to give easy access for downstream processing.

For more detailed instructions, run :code:`multiqc -h` or see the MultiQC website at http://multiqc.info
"""

from setuptools import setup, find_packages

version = '0.3.2dev'

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
        'pyyaml',
        'click'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'qualimap = multiqc.modules.qualimap:MultiqcModule',
            'featureCounts = multiqc.modules.featureCounts:MultiqcModule',
            'preseq = multiqc.modules.preseq:MultiqcModule',
            'picard = multiqc.modules.picard:MultiqcModule',
            'bismark = multiqc.modules.bismark:MultiqcModule',
            'star = multiqc.modules.star:MultiqcModule',
            'tophat = multiqc.modules.tophat:MultiqcModule',
            'bowtie2 = multiqc.modules.bowtie2:MultiqcModule',
            'bowtie1 = multiqc.modules.bowtie1:MultiqcModule',
            'cutadapt = multiqc.modules.cutadapt:MultiqcModule',
            'fastq_screen = multiqc.modules.fastq_screen:MultiqcModule',
            'fastqc = multiqc.modules.fastqc:MultiqcModule',
        ],
        'multiqc.templates.v1': [
            'default = multiqc.templates.default',
            'default_dev = multiqc.templates.default_dev',
            'geo = multiqc.templates.geo',
        ],
        # 'multiqc.hooks.v1': [
            # 'execution_start = myplugin.hooks:execution_start',
            # 'before_modules = myplugin.hooks:before_modules',
            # 'after_modules = myplugin.hooks:after_modules',
            # 'execution_finish = myplugin.hooks:execution_finish',
        # ]
    },
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

