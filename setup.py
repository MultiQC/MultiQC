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

MultiQC was written by Phil Ewels (http://phil.ewels.co.uk) at Seqera Labs (https://seqera.io/), originally at SciLifeLab Sweden (http://www.scilifelab.se)
"""

from setuptools import find_packages, setup

version = "1.18dev"
dl_version = "master" if "dev" in version else "v{}".format(version)

print(
    f"""-----------------------------------
 Installing MultiQC version {version}
-----------------------------------

"""
)

setup(
    name="multiqc",
    version=version,
    author="Phil Ewels",
    author_email="phil.ewels@seqera.io",
    description="Create aggregate bioinformatics analysis reports across many samples and tools",
    long_description=__doc__,
    keywords=["bioinformatics", "biology", "sequencing", "NGS", "next generation sequencing", "quality control"],
    url="http://multiqc.info",
    download_url="https://github.com/ewels/MultiQC/tarball/{}".format(dl_version),
    license="GPLv3",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        "matplotlib>=2.1.1",
        "networkx>=2.5.1",
        "numpy",
        "click",
        "coloredlogs",
        "future>0.14.0",
        "jinja2>=3.0.0",
        "markdown",
        "packaging",
        "pyyaml>=4",
        "requests",
        "rich>=10",
        "rich-click",
        "spectra>=0.0.10",
        "importlib-metadata",
    ],
    entry_points={
        "console_scripts": [
            "multiqc=multiqc.__main__:run_multiqc",
        ],
        "multiqc.modules.v1": [
            "adapterRemoval = multiqc.modules.adapterRemoval.adapterRemoval:MultiqcModule",
            "afterqc = multiqc.modules.afterqc.afterqc:MultiqcModule",
            "anglerfish = multiqc.modules.anglerfish.anglerfish:MultiqcModule",
            "bakta = multiqc.modules.bakta.bakta:MultiqcModule",
            "bamtools = multiqc.modules.bamtools.bamtools:MultiqcModule",
            "bbduk = multiqc.modules.bbduk.bbduk:MultiqcModule",
            "bbmap = multiqc.modules.bbmap.bbmap:MultiqcModule",
            "bcftools = multiqc.modules.bcftools.bcftools:MultiqcModule",
            "bcl2fastq = multiqc.modules.bcl2fastq.bcl2fastq:MultiqcModule",
            "bclconvert = multiqc.modules.bclconvert.bclconvert:MultiqcModule",
            "biobambam2 = multiqc.modules.biobambam2.biobambam2:MultiqcModule",
            "biobloomtools = multiqc.modules.biobloomtools.biobloomtools:MultiqcModule",
            "biscuit = multiqc.modules.biscuit.biscuit:MultiqcModule",
            "bismark = multiqc.modules.bismark.bismark:MultiqcModule",
            "bowtie1 = multiqc.modules.bowtie1.bowtie1:MultiqcModule",
            "bowtie2 = multiqc.modules.bowtie2.bowtie2:MultiqcModule",
            "bracken = multiqc.modules.bracken.bracken:MultiqcModule",
            "busco = multiqc.modules.busco.busco:MultiqcModule",
            "bustools = multiqc.modules.bustools.bustools:MultiqcModule",
            "ccs = multiqc.modules.ccs:.:MultiqcModule",
            "cellranger = multiqc.modules.cellranger.cellranger:MultiqcModule",
            "checkqc = multiqc.modules.checkqc.checkqc:MultiqcModule",
            "clipandmerge = multiqc.modules.clipandmerge.clipandmerge:MultiqcModule",
            "clusterflow = multiqc.modules.clusterflow.clusterflow:MultiqcModule",
            "conpair = multiqc.modules.conpair.conpair:MultiqcModule",
            "custom_content = multiqc.modules.custom_content.custom_content:custom_module_classes",  # special case
            "cutadapt = multiqc.modules.cutadapt.cutadapt:MultiqcModule",
            "damageprofiler = multiqc.modules.damageprofiler.damageprofiler:MultiqcModule",
            "dedup = multiqc.modules.dedup.dedup:MultiqcModule",
            "deeptools = multiqc.modules.deeptools.deeptools:MultiqcModule",
            "diamond = multiqc.modules.diamond.diamond:MultiqcModule",
            "disambiguate = multiqc.modules.disambiguate.disambiguate:MultiqcModule",
            "dragen = multiqc.modules.dragen.dragen:MultiqcModule",
            "dragen_fastqc = multiqc.modules.dragen_fastqc.dragen_fastqc:MultiqcModule",
            "eigenstratdatabasetools = multiqc.modules.eigenstratdatabasetools.eigenstratdatabasetools:MultiqcModule",
            "fastp = multiqc.modules.fastp.fastp:MultiqcModule",
            "fastq_screen = multiqc.modules.fastq_screen.fastq_screen:MultiqcModule",
            "fastqc = multiqc.modules.fastqc.fastqc:MultiqcModule",
            "featureCounts = multiqc.modules.featureCounts.featureCounts:MultiqcModule",
            "fgbio = multiqc.modules.fgbio.fgbio:MultiqcModule",
            "flash = multiqc.modules.flash.flash:MultiqcModule",
            "flexbar = multiqc.modules.flexbar.flexbar:MultiqcModule",
            "filtlong = multiqc.modules.filtlong.filtlong:MultiqcModule",
            "freyja = multiqc.modules.freyja.freyja:MultiqcModule",
            "gffcompare = multiqc.modules.gffcompare.gffcompare:MultiqcModule",
            "gatk = multiqc.modules.gatk:.:MultiqcModule",
            "goleft_indexcov = multiqc.modules.goleft_indexcov.goleft_indexcov:MultiqcModule",
            "gopeaks = multiqc.modules.gopeaks.gopeaks:MultiqcModule",
            "happy = multiqc.modules.happy.happy:MultiqcModule",
            "hicexplorer = multiqc.modules.hicexplorer.hicexplorer:MultiqcModule",
            "hicpro = multiqc.modules.hicpro.hicpro:MultiqcModule",
            "hicup = multiqc.modules.hicup.hicup:MultiqcModule",
            "hifiasm = multiqc.modules.hifiasm.hifiasm:MultiqcModule",
            "hisat2 = multiqc.modules.hisat2.hisat2:MultiqcModule",
            "homer = multiqc.modules.homer.homer:MultiqcModule",
            "hops = multiqc.modules.hops:.:MultiqcModule",
            "htseq = multiqc.modules.htseq.htseq:MultiqcModule",
            "humid = multiqc.modules.humid.humid:MultiqcModule",
            "interop = multiqc.modules.interop.interop:MultiqcModule",
            "ivar = multiqc.modules.ivar:.:MultiqcModule",
            "jcvi = multiqc.modules.jcvi:.:MultiqcModule",
            "jellyfish = multiqc.modules.jellyfish.jellyfish:MultiqcModule",
            "kaiju = multiqc.modules.kaiju.kaiju:MultiqcModule",
            "kallisto = multiqc.modules.kallisto.kallisto:MultiqcModule",
            "kat = multiqc.modules.kat:.:MultiqcModule",
            "kraken = multiqc.modules.kraken.kraken:MultiqcModule",
            "leehom = multiqc.modules.leehom.leehom:MultiqcModule",
            "librarian = multiqc.modules.librarian.librarian:MultiqcModule",
            "lima = multiqc.modules.lima:.:MultiqcModule",
            "longranger = multiqc.modules.longranger.longranger:MultiqcModule",
            "macs2 = multiqc.modules.macs2.macs2:MultiqcModule",
            "malt = multiqc.modules.malt:.:MultiqcModule",
            "mapdamage = multiqc.modules.mapdamage.mapdamage:MultiqcModule",
            "methylQA = multiqc.modules.methylQA.methylQA:MultiqcModule",
            "minionqc = multiqc.modules.minionqc.minionqc:MultiqcModule",
            "mirtop = multiqc.modules.mirtop.mirtop:MultiqcModule",
            "mirtrace = multiqc.modules.mirtrace.mirtrace:MultiqcModule",
            "mosdepth = multiqc.modules.mosdepth.mosdepth:MultiqcModule",
            "motus = multiqc.modules.motus.motus:MultiqcModule",
            "mtnucratio = multiqc.modules.mtnucratio.mtnucratio:MultiqcModule",
            "multivcfanalyzer = multiqc.modules.multivcfanalyzer.multivcfanalyzer:MultiqcModule",
            "nanostat = multiqc.modules.nanostat.nanostat:MultiqcModule",
            "nextclade = multiqc.modules.nextclade.nextclade:MultiqcModule",
            "ngsderive = multiqc.modules.ngsderive.ngsderive:MultiqcModule",
            "odgi = multiqc.modules.odgi:.:MultiqcModule",
            "optitype = multiqc.modules.optitype.optitype:MultiqcModule",
            "pangolin = multiqc.modules.pangolin.pangolin:MultiqcModule",
            "pbmarkdup = multiqc.modules.pbmarkdup.pbmarkdup:MultiqcModule",
            "peddy = multiqc.modules.peddy.peddy:MultiqcModule",
            "phantompeakqualtools = multiqc.modules.phantompeakqualtools.phantompeakqualtools:MultiqcModule",
            "picard = multiqc.modules.picard.picard:MultiqcModule",
            "porechop = multiqc.modules.porechop.porechop:MultiqcModule",
            "preseq = multiqc.modules.preseq.preseq:MultiqcModule",
            "prinseqplusplus = multiqc.modules.prinseqplusplus.prinseqplusplus:MultiqcModule",
            "prokka = multiqc.modules.prokka.prokka:MultiqcModule",
            "purple = multiqc.modules.purple.purple:MultiqcModule",
            "pychopper = multiqc.modules.pychopper.pychopper:MultiqcModule",
            "pycoqc = multiqc.modules.pycoqc.pycoqc:MultiqcModule",
            "qc3C = multiqc.modules.qc3C:.:MultiqcModule",
            "qorts = multiqc.modules.qorts.qorts:MultiqcModule",
            "qualimap = multiqc.modules.qualimap.qualimap:MultiqcModule",
            "quast = multiqc.modules.quast.quast:MultiqcModule",
            "rna_seqc = multiqc.modules.rna_seqc.rna_seqc:MultiqcModule",
            "rockhopper = multiqc.modules.rockhopper.rockhopper:MultiqcModule",
            "rsem = multiqc.modules.rsem:.:MultiqcModule",
            "rseqc = multiqc.modules.rseqc.rseqc:MultiqcModule",
            "salmon = multiqc.modules.salmon.salmon:MultiqcModule",
            "sambamba = multiqc.modules.sambamba.sambamba:MultiqcModule",
            "samblaster = multiqc.modules.samblaster.samblaster:MultiqcModule",
            "samtools = multiqc.modules.samtools.samtools:MultiqcModule",
            "sargasso = multiqc.modules.sargasso.sargasso:MultiqcModule",
            "sentieon = multiqc.modules.sentieon.sentieon:MultiqcModule",
            "seqyclean = multiqc.modules.seqyclean.seqyclean:MultiqcModule",
            "sexdeterrmine = multiqc.modules.sexdeterrmine.sexdeterrmine:MultiqcModule",
            "sickle = multiqc.modules.sickle.sickle:MultiqcModule",
            "skewer = multiqc.modules.skewer.skewer:MultiqcModule",
            "slamdunk = multiqc.modules.slamdunk.slamdunk:MultiqcModule",
            "snippy = multiqc.modules.snippy.snippy:MultiqcModule",
            "snpeff = multiqc.modules.snpeff.snpeff:MultiqcModule",
            "snpsplit = multiqc.modules.snpsplit.snpsplit:MultiqcModule",
            "somalier = multiqc.modules.somalier.somalier:MultiqcModule",
            "sortmerna = multiqc.modules.sortmerna.sortmerna:MultiqcModule",
            "sourmash = multiqc.modules.sourmash.sourmash:MultiqcModule",
            "stacks = multiqc.modules.stacks.stacks:MultiqcModule",
            "star = multiqc.modules.star:.:MultiqcModule",
            "supernova = multiqc.modules.supernova.supernova:MultiqcModule",
            "theta2 = multiqc.modules.theta2.theta2:MultiqcModule",
            "tophat = multiqc.modules.tophat.tophat:MultiqcModule",
            "trimmomatic = multiqc.modules.trimmomatic.trimmomatic:MultiqcModule",
            "truvari = multiqc.modules.truvari.truvari:MultiqcModule",
            "umitools = multiqc.modules.umitools.umitools:MultiqcModule",
            "varscan2 = multiqc.modules.varscan2.varscan2:MultiqcModule",
            "vcftools = multiqc.modules.vcftools.vcftools:MultiqcModule",
            "vep = multiqc.modules.vep:.:MultiqcModule",
            "verifybamid = multiqc.modules.verifybamid.verifybamid:MultiqcModule",
            "whatshap = multiqc.modules.whatshap.whatshap:MultiqcModule",
        ],
        "multiqc.templates.v1": [
            "default = multiqc.templates.default",
            "default_dev = multiqc.templates.default_dev",
            "sections = multiqc.templates.sections",
            "simple = multiqc.templates.simple",
            "gathered = multiqc.templates.gathered",
            "geo = multiqc.templates.geo",
        ],
        ## See https://multiqc.info/docs/#multiqc-plugins for documentation
        # 'multiqc.cli_options.v1': [
        #     'my-new-option = myplugin.cli:new_option'
        # ],
        # 'multiqc.hooks.v1': [
        #     'before_config = myplugin.hooks:before_config',
        #     'config_loaded = myplugin.hooks:config_loaded',
        #     'execution_start = myplugin.hooks:execution_start',
        #     'before_modules = myplugin.hooks:before_modules',
        #     'after_modules = myplugin.hooks:after_modules',
        #     'before_report_generation = myplugin.hooks:before_report_generation',
        #     'execution_finish = myplugin.hooks:execution_finish',
        # ]
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: JavaScript",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)

print(
    """
--------------------------------
 MultiQC installation complete!
--------------------------------
For help in running MultiQC, please see the documentation available
at http://multiqc.info or run: multiqc --help
"""
)
