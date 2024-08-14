---
name: PURPLE
urls: ["https://github.com/hartwigmedical/hmftools/"]
summary: >
  A purity, ploidy and copy number estimator for whole genome tumor data
---

PURPLE combines B-allele frequency (BAF), read depth ratios, somatic variants and
structural variant breakpoints to estimate the purity and copy number profile
of a tumor sample, and also predicts gender, the MSI status, tumor mutational
load and burden, clonality and the whole genome duplication status.

### File search patterns

```yaml
purple/purity:
  fn: "*.purple.purity.tsv"
purple/qc:
  fn: "*.purple.qc"
```
