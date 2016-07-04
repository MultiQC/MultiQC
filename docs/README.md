---
Using MultiQC:
  Installation: installation.md
  Running MultiQC: usage.md
  Using Reports: reports.md
  Configuration: config.md
  Common Problems: troubleshooting.md
MultiQC Modules:
  Pre-alignment:
    FastQC: modules/fastqc.md
    FastQ Screen: modules/fastq_screen.md
    Skewer: modules/skewer.md
    Cutadapt: modules/cutadapt.md
    Trimmomatic: modules/trimmomatic.md
  Aligners:
    Bowtie 1: modules/bowtie1.md
    Bowtie 2: modules/bowtie2.md
    TopHat: modules/tophat.md
    STAR: modules/star.md
    Salmon: modules/salmon.md
    Kallisto: modules/kallisto.md
    HiCUP: modules/hicup.md
    Bismark: modules/bismark.md
  Post-alignment:
    Bamtools: modules/bamtools.md
    Samtools: modules/samtools.md
    Samblaster: modules/samblaster.md
    Preseq: modules/preseq.md
    Picard: modules/picard.md
    RSeQC: modules/rseqc.md
    Methyl QA: modules/methylQA.md
    featureCounts: modules/featureCounts.md
    Qualimap: modules/qualimap.md
    SnpEff: modules/snpeff.md
Coding with MultiQC:
  Writing new templates: templates.md
  Writing new modules: modules.md
  MultiQC Plugins: plugins.md
  JavaScript Functions: javascript.md
---

# Welcome!

## MultiQC v0.8dev Documentation

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
   - [Troubleshooting](troubleshooting.md)
 - MultiQC Modules
   - [FastQC](fastqc.md)
 - Coding with MultiQC
   - [Writing new templates](templates.md)
   - [Writing new modules](modules.md)
   - [Plugins](plugins.md)
   - [JavaScript](javascript.md)

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
[contributing instructions](https://github.com/ewels/MultiQC/blob/master/CONTRIBUTING.md).
