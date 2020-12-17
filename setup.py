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
import sys

version = '1.10dev'
dl_version = 'master' if 'dev' in version else 'v{}'.format(version)

print("""-----------------------------------
 Installing MultiQC version {}
-----------------------------------

""".format(version))

# Set version requirements according to what version of Python we're running
networkx_version = ''
numpy_version = ''
matplotlib_version = '>=2.1.1'
jinja2_version = '>=2.9'
markdown_version = ''
if sys.version_info[0:2] < (3, 6):
    # Lots of tools have dropped Python 3 support, so limit their versions
    matplotlib_version += ',<3.0.0'
    numpy_version = '<1.17'
    networkx_version = '<2.3'
    jinja2_version += ',<3.0'
    markdown_version = '<3.2'

install_requires = [
    'matplotlib' + matplotlib_version,
    'networkx' + networkx_version,
    'numpy' + numpy_version,
    'click',
    'coloredlogs',
    'future>0.14.0',
    'jinja2' + jinja2_version,
    'lzstring',
    'markdown' + markdown_version,
    'pyyaml>=4',
    'requests',
    'simplejson',
    'spectra>=0.0.10',
]

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
    install_requires = install_requires,
    entry_points = {
        "console_scripts": [
            "multiqc=multiqc.__main__:multiqc",
        ],
        'multiqc.modules.v1': [
            'adapterRemoval = multiqc.modules.adapterRemoval:MultiqcModule',
            'afterqc = multiqc.modules.afterqc:MultiqcModule',
            'artic = multiqc.modules.artic:MultiqcModule',
            'bamtools = multiqc.modules.bamtools:MultiqcModule',
            'bbmap = multiqc.modules.bbmap:MultiqcModule',
            'bcftools = multiqc.modules.bcftools:MultiqcModule',
            'bcl2fastq = multiqc.modules.bcl2fastq:MultiqcModule',
            'biobambam2 = multiqc.modules.biobambam2:MultiqcModule',
            'biobloomtools = multiqc.modules.biobloomtools:MultiqcModule',
            'biscuit = multiqc.modules.biscuit:MultiqcModule',
            'bismark = multiqc.modules.bismark:MultiqcModule',
            'bowtie1 = multiqc.modules.bowtie1:MultiqcModule',
            'bowtie2 = multiqc.modules.bowtie2:MultiqcModule',
            'busco = multiqc.modules.busco:MultiqcModule',
            'clipandmerge = multiqc.modules.clipandmerge:MultiqcModule',
            'clusterflow = multiqc.modules.clusterflow:MultiqcModule',
            'conpair = multiqc.modules.conpair:MultiqcModule',
            'custom_content = multiqc.modules.custom_content:custom_module_classes', # special case
            'cutadapt = multiqc.modules.cutadapt:MultiqcModule',
            'damageprofiler = multiqc.modules.damageprofiler:MultiqcModule',
            'dedup = multiqc.modules.dedup:MultiqcModule',
            'deeptools = multiqc.modules.deeptools:MultiqcModule',
            'disambiguate = multiqc.modules.disambiguate:MultiqcModule',
            'dragen = multiqc.modules.dragen:MultiqcModule',
            'fastp = multiqc.modules.fastp:MultiqcModule',
            'fastq_screen = multiqc.modules.fastq_screen:MultiqcModule',
            'fastqc = multiqc.modules.fastqc:MultiqcModule',
            'featureCounts = multiqc.modules.featureCounts:MultiqcModule',
            'fgbio = multiqc.modules.fgbio:MultiqcModule',
            'flash = multiqc.modules.flash:MultiqcModule',
            'flexbar = multiqc.modules.flexbar:MultiqcModule',
            'gatk = multiqc.modules.gatk:MultiqcModule',
            'goleft_indexcov = multiqc.modules.goleft_indexcov:MultiqcModule',
            'happy = multiqc.modules.happy:MultiqcModule',
            'hicexplorer = multiqc.modules.hicexplorer:MultiqcModule',
            'hicpro = multiqc.modules.hicpro:MultiqcModule',
            'hicup = multiqc.modules.hicup:MultiqcModule',
            'hisat2 = multiqc.modules.hisat2:MultiqcModule',
            'homer = multiqc.modules.homer:MultiqcModule',
            'hops = multiqc.modules.hops:MultiqcModule',
            'htseq = multiqc.modules.htseq:MultiqcModule',
            'interop = multiqc.modules.interop:MultiqcModule',
            'ivar = multiqc.modules.ivar:MultiqcModule',
            'jellyfish = multiqc.modules.jellyfish:MultiqcModule',
            'kaiju = multiqc.modules.kaiju:MultiqcModule',
            'kallisto = multiqc.modules.kallisto:MultiqcModule',
            'kat = multiqc.modules.kat:MultiqcModule',
            'kraken = multiqc.modules.kraken:MultiqcModule',
            'leehom = multiqc.modules.leehom:MultiqcModule',
            'longranger = multiqc.modules.longranger:MultiqcModule',
            'macs2 = multiqc.modules.macs2:MultiqcModule',
            'malt = multiqc.modules.malt:MultiqcModule',
            'methylQA = multiqc.modules.methylQA:MultiqcModule',
            'minionqc = multiqc.modules.minionqc:MultiqcModule',
            'mirtop = multiqc.modules.mirtop:MultiqcModule',
            'mirtrace = multiqc.modules.mirtrace:MultiqcModule',
            'mosdepth = multiqc.modules.mosdepth:MultiqcModule',
            'mtnucratio = multiqc.modules.mtnucratio:MultiqcModule',
            'multivcfanalyzer = multiqc.modules.multivcfanalyzer:MultiqcModule',
            'peddy = multiqc.modules.peddy:MultiqcModule',
            'phantompeakqualtools = multiqc.modules.phantompeakqualtools:MultiqcModule',
            'picard = multiqc.modules.picard:MultiqcModule',
            'preseq = multiqc.modules.preseq:MultiqcModule',
            'prokka = multiqc.modules.prokka:MultiqcModule',
            'pychopper = multiqc.modules.pychopper:MultiqcModule',
            'purple = multiqc.modules.purple:MultiqcModule',
            'pycoqc = multiqc.modules.pycoqc:MultiqcModule',
            'qc3C = multiqc.modules.qc3C:MultiqcModule',
            'qorts = multiqc.modules.qorts:MultiqcModule',
            'qualimap = multiqc.modules.qualimap:MultiqcModule',
            'quast = multiqc.modules.quast:MultiqcModule',
            'rna_seqc = multiqc.modules.rna_seqc:MultiqcModule',
            'rockhopper = multiqc.modules.rockhopper:MultiqcModule',
            'rsem = multiqc.modules.rsem:MultiqcModule',
            'rseqc = multiqc.modules.rseqc:MultiqcModule',
            'salmon = multiqc.modules.salmon:MultiqcModule',
            'samblaster = multiqc.modules.samblaster:MultiqcModule',
            'samtools = multiqc.modules.samtools:MultiqcModule',
            'sargasso = multiqc.modules.sargasso:MultiqcModule',
            'sentieon = multiqc.modules.sentieon:MultiqcModule',
            'seqyclean = multiqc.modules.seqyclean:MultiqcModule',
            'sexdeterrmine = multiqc.modules.sexdeterrmine:MultiqcModule',
            'sickle = multiqc.modules.sickle:MultiqcModule',
            'skewer = multiqc.modules.skewer:MultiqcModule',
            'slamdunk = multiqc.modules.slamdunk:MultiqcModule',
            'snpeff = multiqc.modules.snpeff:MultiqcModule',
            'snpsplit = multiqc.modules.snpsplit:MultiqcModule',
            'somalier = multiqc.modules.somalier:MultiqcModule',
            'sortmerna = multiqc.modules.sortmerna:MultiqcModule',
            'stacks = multiqc.modules.stacks:MultiqcModule',
            'star = multiqc.modules.star:MultiqcModule',
            'supernova = multiqc.modules.supernova:MultiqcModule',
            'theta2 = multiqc.modules.theta2:MultiqcModule',
            'tophat = multiqc.modules.tophat:MultiqcModule',
            'trimmomatic = multiqc.modules.trimmomatic:MultiqcModule',
            'varscan2 = multiqc.modules.varscan2:MultiqcModule',
            'vcftools = multiqc.modules.vcftools:MultiqcModule',
            'verifybamid = multiqc.modules.verifybamid:MultiqcModule'
        ],
        'multiqc.templates.v1': [
            'default = multiqc.templates.default',
            'default_dev = multiqc.templates.default_dev',
            'sections = multiqc.templates.sections',
            'simple = multiqc.templates.simple',
            'gathered = multiqc.templates.gathered',
            'geo = multiqc.templates.geo',
        ],
## See https://multiqc.info/docs/#multiqc-plugins for documentation
#       'multiqc.cli_options.v1': [
#           'my-new-option = myplugin.cli:new_option'
#       ],
#       'multiqc.hooks.v1': [
#           'before_config = myplugin.hooks:before_config',
#           'config_loaded = myplugin.hooks:config_loaded',
#           'execution_start = myplugin.hooks:execution_start',
#           'before_modules = myplugin.hooks:before_modules',
#           'after_modules = myplugin.hooks:after_modules',
#           'before_report_generation = myplugin.hooks:before_report_generation',
#           'execution_finish = myplugin.hooks:execution_finish',
#       ]
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
