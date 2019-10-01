---
Name: ABRicate
URL: https://github.com/tseemann/abricate
Description: >
    Mass screening of contigs for antimicrobial resistance or virulence genes
---

ABRicate uses a blast search against a prepared database to locate genes in fasta files.

This modules parses the tab-delimited summary output of abricate into a heatmap showing the presence and absence of genes identified.

Example usage: `abricate --summary 1.tab 2.tab > abricate_summary.txt`
