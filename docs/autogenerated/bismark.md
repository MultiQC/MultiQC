---
name: Bismark
urls: ['http://www.bioinformatics.babraham.ac.uk/projects/bismark/']
summary: >
  Maps bisulfite converted sequence reads and determine cytosine methylation states
---

### File search patterns

```yaml
bismark/align:
  fn: '*_[SP]E_report.txt'
bismark/bam2nuc:
  fn: '*.nucleotide_stats.txt'
bismark/dedup:
  fn: '*.deduplication_report.txt'
bismark/m_bias:
  fn: '*M-bias.txt'
bismark/meth_extract:
  fn: '*_splitting_report.txt'
```