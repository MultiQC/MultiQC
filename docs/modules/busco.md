---
Name: BUSCO
URL: http://busco.ezlab.org/
Description: >
    BUSCO assesses genome assembly and annotation completeness with
    Benchmarking Universal Single-Copy Orthologs.
---

BUSCO v2 provides quantitative measures for the assessment of genome
assembly, gene set, and transcriptome completeness, based on
evolutionarily-informed expectations of gene content from near-universal
single-copy orthologs selected from OrthoDB v9.

The MultiQC module parses the `short_summary_[samplename].txt` files and
plots the proportion of BUSCO types found. MultiQC has been tested with
output from BUSCO v1.22 - v2.
