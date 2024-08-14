---
name: HiC-Pro
urls: ["https://github.com/nservant/HiC-Pro"]
summary: >
  Pipeline for Hi-C data processing
---

**Note** - because this module shares sample identifiers across multiple files,
the `--fn_as_s_name` / `config.use_filename_as_sample_name` functionality has been disabled and has no effect.

The MultiQC module is supported since HiC-Pro v2.11.0.

### File search patterns

```yaml
hicpro/assplit:
  fn: "*assplit.stat"
hicpro/mRSstat:
  contents: Valid_interaction_pairs
  fn: "*RSstat"
hicpro/mergestat:
  contents: valid_interaction
  fn: "*.mergestat"
  num_lines: 10
hicpro/mmapstat:
  contents: total_R
  fn: "*mapstat"
  num_lines: 10
hicpro/mpairstat:
  contents: Total_pairs_processed
  fn: "*pairstat"
  num_lines: 10
```
