---
name: miRTrace
urls: ["https://github.com/friedlanderlab/mirtrace"]
summary: >
  Quality control for small RNA sequencing data
---

miRTrace performs adapter trimming and discards the reads that fail to pass
the QC filters. miRTrace specifically addresses sequencing quality, read length,
sequencing depth and miRNA complexity and also identifies the presence of both
miRNAs and undesirable sequences derived from tRNAs, rRNAs, or Illumina artifact
sequences.
miRTrace also profiles clade-specific miRNAs based on a comprehensive catalog
of clade-specific miRNA families identified previously. With this information,
miRTrace can detect exogenous miRNAs, which could be contamination derived,
e.g. index mis-assignment on sample demultiplexing, or biologically derived,
e.g. parasitic RNAs.

### File search patterns

```yaml
mirtrace/contaminationbasic:
  fn: mirtrace-stats-contamination_basic.tsv
mirtrace/length:
  fn: mirtrace-stats-length.tsv
mirtrace/mirnacomplexity:
  fn: mirtrace-stats-mirna-complexity.tsv
mirtrace/summary:
  fn: mirtrace-results.json
```