---
Name: Sambamba
URL: https://lomereiter.github.io/sambamba/
Description: >
  Sambamba is a suite of programs written in the D Language for users to 
  process high-throughput sequencing data.
---

[Sambamba](https://lomereiter.github.io/sambamba/) is a suite of programs for
users to quickly and efficiently process their high-throughput sequencing data.
It is functionally similar to Samtools, but the source code is written in the
D Language; it allows for faster performance while still being easy to use.

Supported commands:

- `markdup`

### markdup

This module parses key phrases in the output log files to find duplicate +
unique reads and then calculates duplicate rate per sample. It will
will work for both single and paired-end data.
The absolute number of reads by type are displayed in a stacked bar plot,
and duplicate rates are in the general statistics table.

Duplicate rates are calculated as follows:

#### Paired end

> `duplicate_rate = duplicateReads / (sortedEndPairs * 2 + singleEnds - singleUnmatchedPairs) * 100`

#### Single end

> `duplicate_rate = duplicateReads / singleEnds * 100`

If Sambamba Markdup is invoked using Snakemake, the following bare-bones
rule should work fine:

```python
rule markdup:
  input:
    "data/align/{sample}.bam"
  output:
    "data/markdup/{sample}.markdup.bam"
  log:
    "data/logs/{sample}.log"
  shell:
    "sambamba markdup {input} {output} > {log} 2>&1"
```
