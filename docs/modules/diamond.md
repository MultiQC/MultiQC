---
name: DIAMOND
url: https://github.com/bbuchfink/diamond
description: >
  DIAMOND is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data.
---

Sequence aligner for protein and translated DNA searches and functions as a drop-in replacement for the NCBI BLAST software tools. The key features are:

- Pairwise alignment of proteins and translated DNA at 100x-10,000x speed of BLAST.
- Frameshift alignments for long read analysis.
- Low resource requirements and suitable for running on standard desktops or laptops.
- Various output formats, including BLAST pairwise, tabular and XML, as well as taxonomic classification.

The module takes summary statistics from the `diamond.log` file (`--log` option). It parses and reports the number of sequences aligned and displays them in the General Stats table.
