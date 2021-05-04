---
Name: Sambamba
URL: https://lomereiter.github.io/sambamba/
Description: >
  Sambamba is a suite of programs written in the D Language for users to 
  process high-throughput sequencing data.
---

[Sambamba](https://lomereiter.github.io/sambamba/) is a suite of programs for users to quickly and efficiently process their high-throughput sequencing data. It is functionally similar to Samtools, but the source code is written in the D Language; it allows for faster performance while still being easy to use.

Supported commands:

- `markdup`

### markdup

This module parses key phrases in the output log files to find duplicate and non-duplicate reads. The absolute read counts are displayed in a stacked bar plot, and the duplicate rate per sample is calculated using the following formula:

>duplicate_rate = duplicate_reads / (sorted_end_pairs x2 + single_ends - single_unmatched_pairs) x100

The duplicate rates are displayed in the general statistics table at the top of the report. Lastly, If Sambamba Markdup is invoked using Snakemake, the following bare-bones shell command should work fine:

```bash
sambamba markdup {input} {output} > {log} 2>&1
```
