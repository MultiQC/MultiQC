---
name: HOMER
url: http://homer.ucsd.edu/homer/
description: Motif discovery and next-gen sequencing analysis
---

HOMER contains many useful tools for analyzing ChIP-Seq, GRO-Seq, RNA-Seq, DNase-Seq, Hi-C and numerous
other types of functional genomics sequencing data sets. The module currently only parses output from the
`findPeaks` and `TagDirectory` tools. If you would like support to be added for other HOMER tools please
open a [new issue](https://github.com/MultiQC/MultiQC/issues/new) on the MultiQC GitHub page.

#### FindPeaks

The HOMER findPeaks MultiQC module parses the summary statistics found at the top
of HOMER peak files. Three key statistics are shown in the General Statistics table,
all others are saved to `multiqc_data/multiqc_homer_findpeaks.txt`.

#### TagDirectory

The HOMER tag directory submodule parses output from files
[tag directory](http://homer.ucsd.edu/homer/ngs/tagDir.html) output files, generating
a number of diagnostic plots.
