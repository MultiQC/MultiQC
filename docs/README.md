---
Using MultiQC:
  Installation: installation.md
  Running MultiQC: usage.md
  Using Reports: reports.md
  Configuration: config.md
  Customising Reports: customisation.md
  Common Problems: troubleshooting.md
MultiQC Modules:
  Pre-alignment:
    Adapter Removal: modules/adapterRemoval.md
    Cluster Flow: modules/clusterflow.md
    Cutadapt: modules/cutadapt.md
    FastQC: modules/fastqc.md
    FastQ Screen: modules/fastq_screen.md
    Skewer: modules/skewer.md
    SortMeRNA: modules/sortmerna.md
    Trimmomatic: modules/trimmomatic.md
  Aligners:
    Bismark: modules/bismark.md
    Bowtie 1: modules/bowtie1.md
    Bowtie 2: modules/bowtie2.md
    HiCUP: modules/hicup.md
    Kallisto: modules/kallisto.md
    STAR: modules/star.md
    Salmon: modules/salmon.md
    TopHat: modules/tophat.md
  Post-alignment:
    Bamtools: modules/bamtools.md
    Bcftools: modules/bcftools.md
    BUSCO: modules/busco.md
    featureCounts: modules/featureCounts.md
    GATK: modules/gatk.md
    goleft_indexcov: modules/goleft_indexcov.md
    HTSeq: modules/htseq.md
    Methyl QA: modules/methylQA.md
    Peddy: modules/peddy.md
    Picard: modules/picard.md
    Preseq: modules/preseq.md
    Prokka: modules/prokka.md
    Qualimap: modules/qualimap.md
    Quast: modules/quast.md
    RNA-SeQC: modules/rna_seqc.md
    RSeQC: modules/rseqc.md
    Samblaster: modules/samblaster.md
    Samtools: modules/samtools.md
    Slamdunk: modules/slamdunk.md
    SnpEff: modules/snpeff.md
Custom Content:
  Introduction: custom_content.md
Coding with MultiQC:
  Writing new templates: templates.md
  Writing new modules: modules.md
  Plotting Functions: plots.md
  MultiQC Plugins: plugins.md
  Updating for compatibility: compatibility.md
---

# Welcome!

## MultiQC v1.0dev Documentation

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
