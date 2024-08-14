---
name: PRINSEQ++
urls: ["https://github.com/Adrian-Cantu/PRINSEQ-plus-plus"]
summary: >
  C++ implementation of the prinseq-lite.pl program. Filters, reformats, and trims genomic and metagenomic reads
---

This module requires that PRINSEQ++ has been run with the flag `-VERBOSE 1`.

It uses the log file name as the sample name.

### File search patterns

```yaml
prinseqplusplus:
  - contents: reads removed by -
    num_lines: 2
```
