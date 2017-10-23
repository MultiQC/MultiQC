---
Name: HOMER
URL: http://homer.ucsd.edu/homer/
Description: >
    HOMER is a suite of tools for Motif Discovery and next-gen sequencing analysis.
---

HOMER _(Hypergeometric Optimization of Motif EnRichment)_ is a suite of tools for
Motif Discovery and next-gen sequencing analysis. HOMER contains many useful tools
for analyzing ChIP-Seq, GRO-Seq, RNA-Seq, DNase-Seq, Hi-C and numerous other types
of functional genomics sequencing data sets.

The HOMER MultiQC module currently only parses output from the `findPeaks` tool.
If you would like support to be added for other HOMER tools, please open a
[new issue](https://github.com/ewels/MultiQC/issues/new) on the MultiQC GitHub page.

#### FindPeaks
The HOMER findPeaks MultiQC module parses the summary statistics found at the top
of HOMER peak files. Three key statistics are shown in the General Statistics table,
all others are saved to `multiqc_data/multiqc_homer_findpeaks.txt`.

#### TagDirectory
The HOMER tag directory submodule parses output from files
[tag directory](http://homer.ucsd.edu/homer/ngs/tagDir.html) output files, generating
a number of diagnostic plots.
