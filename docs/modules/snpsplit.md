---
name: SNPsplit
url: https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/
description: SNPsplit is an allele-specific alignment sorter, which is designed to read in alignment files in SAM/BAM format and determine the allelic origin of reads that cover known SNP positions.
---

Currently only the "Allele-tagging" and "Allele-sorting" reports are supported.
The log files from the genome creation steps are not parsed and there are no plots/tables produced from the "SNP coverage" report.

Differences between the numbers in the tagging and sorting reports are due to paired-end reads.
For these, if only a single mate in a pair is assigned to a genome then it will "rescue" its mate and both will be "sorted" into that genome (even though only one of them was tagged).
Conversely, if the mates in a pair are tagged as arising from different genomes, then the pair as a whole is unassignable.
