---
name: HOMER
urls: ['http://homer.ucsd.edu/homer/']
summary: >
  Motif discovery and next-gen sequencing analysis
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

### File search patterns

```yaml
homer/FreqDistribution:
  fn: petag.FreqDistribution_1000.txt
homer/GCcontent:
  fn: tagGCcontent.txt
homer/LengthDistribution:
  fn: tagLengthDistribution.txt
homer/RestrictionDistribution:
  fn: petagRestrictionDistribution.*.txt
homer/findpeaks:
  contents: '# HOMER Peaks'
  num_lines: 3
homer/genomeGCcontent:
  fn: genomeGCcontent.txt
homer/tagInfo:
  fn: tagInfo.txt
```