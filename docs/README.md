---
Using MultiQC:
  Installation: installation.md
  Running MultiQC: usage.md
  Using Reports: reports.md
  Configuration: config.md
  Customising Reports: customisation.md
  Using MultiQC in pipelines: pipelines.md
  Common Problems: troubleshooting.md
MultiQC Modules:
  Pre-alignment:
    Adapter Removal: modules/adapterRemoval.md
    AfterQC: modules/afterqc.md
    Bcl2fastq: modules/bcl2fastq.md
    BioBloom Tools: modules/biobloomtools.md
    Cluster Flow: modules/clusterflow.md
    Cutadapt: modules/cutadapt.md
    ClipAndMerge: modules/clipandmerge.md
    FastQ Screen: modules/fastq_screen.md
    FastQC: modules/fastqc.md
    Fastp: modules/fastp.md
    FLASh: modules/flash.md
    Flexbar: modules/flexbar.md
    InterOp: modules/interop.md
    iVar: modules/ivar.md
    Jellyfish: modules/jellyfish.md
    KAT: modules/kat.md
    leeHom: modules/leehom.md
    minionqc: modules/minionqc.md
    pycoQC: modules/pycoqc.md
    Pychopper: modules/pychopper.md
    SeqyClean: modules/seqyclean.md
    Sickle: modules/sickle.md
    Skewer: modules/skewer.md
    SortMeRNA: modules/sortmerna.md
    Trimmomatic: modules/trimmomatic.md
  Aligners:
    Biscuit: modules/biscuit.md
    Bismark: modules/bismark.md
    Bowtie 1: modules/bowtie1.md
    Bowtie 2: modules/bowtie2.md
    BBMap: modules/bbmap.md
    HiCUP: modules/hicup.md
    HiCPro: modules/hicpro.md
    HISAT2: modules/hisat2.md
    Kallisto: modules/kallisto.md
    Long Ranger: modules/longranger.md
    Salmon: modules/salmon.md
    STAR: modules/star.md
    TopHat: modules/tophat.md
    DRAGEN: modules/dragen.md
    MALT: modules/malt.md
  Post-alignment:
    ARTIC: modules/artic.md
    Bamtools: modules/bamtools.md
    Bcftools: modules/bcftools.md
    biobambam2: modules/biobambam2.md
    BUSCO: modules/busco.md
    Conpair: modules/conpair.md
    DamageProfiler: modules/damageprofiler.md
    DeDup: modules/dedup.md
    deepTools: modules/deeptools.md
    Disambiguate: modules/disambiguate.md
    featureCounts: modules/featureCounts.md
    fgbio: modules/fgbio.md
    GATK: modules/gatk.md
    goleft_indexcov: modules/goleft_indexcov.md
    Hap.py: modules/happy.md
    HiCExplorer: modules/hicexplorer.md
    HOPS: modules/hops.md
    HOMER: modules/homer.md
    HTSeq: modules/htseq.md
    Kaiju: modules/kaiju.md
    Kraken: modules/kraken.md
    MACS2: modules/macs2.md
    Methyl QA: modules/methylQA.md
    mosdepth: modules/mosdepth.md
    miRTrace: modules/mirtrace.md
    mirtop: modules/mirtop.md
    MTNucRatio: modules/mtnucratio.md
    MultiVCFAnalyzer: modules/multivcfanalyzer.md
    phantompeakqualtools: modules/phantompeakqualtools.md
    Peddy: modules/peddy.md
    Picard: modules/picard.md
    Preseq: modules/preseq.md
    Prokka: modules/prokka.md
    PURPLE: modules/purple.md
    qc3C: modules/qc3C.md
    QoRTs: modules/qorts.md
    Qualimap: modules/qualimap.md
    Quast: modules/quast.md
    RNA-SeQC: modules/rna_seqc.md
    Rockhopper: modules/rockhopper.md
    RSEM: modules/rsem.md
    RSeQC: modules/rseqc.md
    Samblaster: modules/samblaster.md
    Samtools: modules/samtools.md
    Sentieon: modules/sentieon.md
    Sargasso: modules/sargasso.md
    SexDetErrmine: modules/sexdeterrmine.md
    Slamdunk: modules/slamdunk.md
    SnpEff: modules/snpeff.md
    SNPsplit: modules/snpsplit.md
    Somalier: modules/somalier.md
    Supernova: modules/supernova.md
    Stacks: modules/stacks.md
    THeTA2: modules/theta2.md
    VarScan2: modules/varscan2.md
    VCFTools: modules/vcftools.md
    verifyBAMID: modules/verifybamid.md
Custom Content:
  Introduction: custom_content.md
Coding with MultiQC:
  Writing new modules: modules.md
  Plotting Functions: plots.md
  MultiQC Plugins: plugins.md
  Writing new templates: templates.md
  Updating for compatibility: compatibility.md
---

# Welcome!

## MultiQC Documentation

MultiQC is a tool to aggregate bioinformatics results across many samples
into a single report. It's written in Python and contains modules for a number
of common tools.

The documentation has the following pages:

 - [Docs homepage](README.md) _(this README file)_
 - Using MultiQC
   - [Installing MultiQC](installation.md)
   - [Running MultiQC](usage.md)
   - [Using Reports](reports.md)
   - [Configuration](config.md)
   - [Customising Reports](customisation.md)
   - [Common Problems](troubleshooting.md)
 - [MultiQC Modules](modules/)
 - [Custom Content](custom_content.md)
 - Coding with MultiQC
   - [Writing new templates](templates.md)
   - [Writing new modules](modules.md)
   - [Plugins](plugins.md)
   - [MultiQC Plugins](plugins.md)
   - [Updating for compatibility](compatibility.md)

These docs can be read in any of three ways:
 - On the MultiQC Website: http://multiqc.info
 - On GitHub: https://github.com/ewels/MultiQC/
 - As part of the distributed source code (in `/docs/`)

If you're curious how the website works, check out the
[MultiQC website repository](https://github.com/ewels/MultiQC_website).

## Contributing to MultiQC

If you write a module which could be of use to others, it would be great to
merge those changes back into the core MultiQC project.

For instructions on how best to do this, please see the
[contributing instructions](https://github.com/ewels/MultiQC/blob/master/.github/CONTRIBUTING.md).
