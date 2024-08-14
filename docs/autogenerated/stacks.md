---
name: Stacks
urls: ["http://catchenlab.life.illinois.edu/stacks/"]
summary: >
  Analyzes restriction enzyme-based data (e.g. RAD-seq)
---

This module is designed to only parse some of the output from the Stacks `denovo_map` pipeline.

The module works with Stacks version 2.1 or greater.

If you are missing some functionality, please submit an issue on the [MultiQC github page](https://github.com/MultiQC/MultiQC)

### File search patterns

```yaml
stacks/gstacks:
  contents: BEGIN effective_coverages_per_sample
  fn: gstacks.log.distribs
stacks/populations:
  contents: BEGIN missing_samples_per_loc_prefilters
  fn: populations.log.distribs
stacks/sumstats:
  contents: "# Pop ID\tPrivate\tNum_Indv\tVar\tStdErr\tP\tVar"
  fn: "*.sumstats_summary.tsv"
  max_filesize: 1000000
```
