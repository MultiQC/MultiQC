#!/usr/bin/env python
"""
MultiQC is a tool to aggregate bioinformatics results across many samples into a single report. It is written in Python and contains modules for a large number of common bioinformatics tools.

You can install MultiQC from PyPI as follows::

    pip install multiqc

Then it's just a case of going to your analysis directory and running the script::

    multiqc .

MultiQC will scan the specified directory (:code:`'.'` is the current dir) and produce a report detailing whatever it finds.

The report is created in :code:`multiqc_report.html` by default. Tab-delimited data files are created in :code:`multiqc_data/` to give easy access for downstream processing.

For more detailed instructions, run :code:`multiqc -h`

See the MultiQC website for documentation and tutorial videos: http://multiqc.info

MultiQC was written by Phil Ewels (http://phil.ewels.co.uk) at SciLifeLab Sweden (http://www.scilifelab.se)
"""

from setuptools import setup, find_packages

version = '1.1'
dl_version = 'master' if 'dev' in version else 'v{}'.format(version)

print("""-----------------------------------
 Installing MultiQC version {}
-----------------------------------

""".format(version))

setup(
    name = 'multiqc',
    version = version,
    author = 'Phil Ewels',
    author_email = 'phil.ewels@scilifelab.se',
    description = "Create aggregate bioinformatics analysis reports across many samples and tools",
    long_description = __doc__,
    keywords = ['bioinformatics', 'biology', 'sequencing', 'NGS', 'next generation sequencing', 'quality control'],
    url = 'http://multiqc.info',
    download_url = 'https://github.com/ewels/MultiQC/tarball/{}'.format(dl_version),
    license = 'GPLv3',
    packages = find_packages(),
    include_package_data = True,
    zip_safe = False,
    scripts = ['scripts/multiqc'],
    install_requires = [
        'click',
        'future>0.14.0',
        'lzstring',
        'jinja2>=2.9',
        'matplotlib',
        'numpy',
        'pyyaml',
        'simplejson',
        'spectra',
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'custom_content = multiqc.modules.custom_content:custom_module_classes', # special case
            'adapterRemoval = multiqc.modules.adapterRemoval:MultiqcModule',
            'bamtools = multiqc.modules.bamtools:MultiqcModule',
            'bcftools = multiqc.modules.bcftools:MultiqcModule',
            'biobloomtools = multiqc.modules.biobloomtools:MultiqcModule',
            'bismark = multiqc.modules.bismark:MultiqcModule',
            'bowtie2 = multiqc.modules.bowtie2:MultiqcModule',
            'bowtie1 = multiqc.modules.bowtie1:MultiqcModule',
            'busco = multiqc.modules.busco:MultiqcModule',
            'clusterflow = multiqc.modules.clusterflow:MultiqcModule',
            'conpair = multiqc.modules.conpair:MultiqcModule',
            'cutadapt = multiqc.modules.cutadapt:MultiqcModule',
            'disambiguate = multiqc.modules.disambiguate:MultiqcModule',
            'fastq_screen = multiqc.modules.fastq_screen:MultiqcModule',
            'fastqc = multiqc.modules.fastqc:MultiqcModule',
            'featureCounts = multiqc.modules.featureCounts:MultiqcModule',
            'flexbar = multiqc.modules.flexbar:MultiqcModule',
            'gatk = multiqc.modules.gatk:MultiqcModule',
            'goleft_indexcov = multiqc.modules.goleft_indexcov:MultiqcModule',
            'hicup = multiqc.modules.hicup:MultiqcModule',
            'hisat2 = multiqc.modules.hisat2:MultiqcModule',
            'homer = multiqc.modules.homer:MultiqcModule',
            'htseq = multiqc.modules.htseq:MultiqcModule',
            'jellyfish = multiqc.modules.jellyfish:MultiqcModule',
            'kallisto = multiqc.modules.kallisto:MultiqcModule',
            'macs2 = multiqc.modules.macs2:MultiqcModule',
            'methylQA = multiqc.modules.methylQA:MultiqcModule',
            'peddy = multiqc.modules.peddy:MultiqcModule',
            'picard = multiqc.modules.picard:MultiqcModule',
            'preseq = multiqc.modules.preseq:MultiqcModule',
            'prokka = multiqc.modules.prokka:MultiqcModule',
            'qorts = multiqc.modules.qorts:MultiqcModule',
            'qualimap = multiqc.modules.qualimap:MultiqcModule',
            'quast = multiqc.modules.quast:MultiqcModule',
            'rna_seqc = multiqc.modules.rna_seqc:MultiqcModule',
            'rseqc = multiqc.modules.rseqc:MultiqcModule',
            'salmon = multiqc.modules.salmon:MultiqcModule',
            'samblaster = multiqc.modules.samblaster:MultiqcModule',
            'samtools = multiqc.modules.samtools:MultiqcModule',
            'skewer = multiqc.modules.skewer:MultiqcModule',
            'slamdunk = multiqc.modules.slamdunk:MultiqcModule',
            'snpeff = multiqc.modules.snpeff:MultiqcModule',
            'star = multiqc.modules.star:MultiqcModule',
            'sortmerna = multiqc.modules.sortmerna:MultiqcModule',
            'theta2 = multiqc.modules.theta2:MultiqcModule',
            'tophat = multiqc.modules.tophat:MultiqcModule',
            'trimmomatic = multiqc.modules.trimmomatic:MultiqcModule',
        ],
        'multiqc.templates.v1': [
            'default = multiqc.templates.default',
            'default_dev = multiqc.templates.default_dev',
            'sections = multiqc.templates.sections',
            'simple = multiqc.templates.simple',
            'geo = multiqc.templates.geo',
        ],
        # 'multiqc.cli_options.v1': [
            # 'my-new-option = myplugin.cli:new_option'
        # ],
        # 'multiqc.hooks.v1': [
            # 'before_config = myplugin.hooks:before_config',
            # 'config_loaded = myplugin.hooks:config_loaded',
            # 'execution_start = myplugin.hooks:execution_start',
            # 'before_modules = myplugin.hooks:before_modules',
            # 'after_modules = myplugin.hooks:after_modules',
            # 'before_report_generation = myplugin.hooks:before_report_generation',
            # 'execution_finish = myplugin.hooks:execution_finish',
        # ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
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

print("""
--------------------------------
 MultiQC installation complete!
--------------------------------
For help in running MultiQC, please see the documentation available
at http://multiqc.info or run: multiqc --help
""")
