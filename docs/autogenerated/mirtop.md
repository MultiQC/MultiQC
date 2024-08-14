---
name: mirtop
urls: ['https://github.com/miRTop/mirtop/']
summary: >
  Annotates miRNAs and isomiRs and compute general statistics in mirGFF3 format
extra_description: >
  This tool is dedicated to the creation and management of miRNA alignment output using the standardized
  GFF3 format (see [miRTop/mirGFF3](https://github.com/miRTop/mirGFF3)).
  A unified miRNA alignment format allows to easily compare the output of different alignment tools.
  Currently, mirtop can convert into mirGFF3 the outputs of commonly used pipelines, such as seqbuster,
  isomiR-SEA, sRNAbench, Prost! as well as BAM files.
---

### File search patterns

```yaml
mirtop:
  fn: '*_mirtop_stats.log'
```