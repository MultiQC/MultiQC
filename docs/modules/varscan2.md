---
name: VarScan2
url: http://dkoboldt.github.io/varscan/
description: Variant detection in massively parallel sequencing data
---

VarScan is a platform-independent mutation caller for targeted, exome, and whole-genome
resequencing data generated on Illumina, SOLiD, Life/PGM, Roche/454, and similar instruments.
The newest version, VarScan 2, is written in Java, so it runs on most operating systems.

VarScan can be used to detect different types of variation:

- Germline variants (SNPs an dindels) in individual samples or pools of samples.
- Multi-sample variants (shared or private) in multi-sample datasets (with mpileup).
- Somatic mutations, LOH events, and germline variants in tumor-normal pairs.
- Somatic copy number alterations (CNAs) in tumor-normal exome data.

The MultiQC module can read output from `mpileup2cns`, `mpileup2snp` and `mpileup2indel` logfiles.
