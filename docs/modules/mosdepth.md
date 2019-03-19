---
Name: mosdepth
URL: https://github.com/brentp/mosdepth
Description: >
    Mosdepth performs fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing.
    It can output:
    - per-base depth
    - mean per-window depth given a window size--as would be used for CNV calling
    - mean per-region given a BED file of regions
    - a distribution of proportion of bases covered at or above a given threshhold for each chromosome and genome-wide.
    - quantized output that merges adjacent bases as long as they fall in the same coverage bins e.g. (10-20)
    - threshold output to indicate how many bases in each region are covered at the given thresholds.
---

[Mosdepth](https://github.com/brentp/mosdepth/) performs fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing.

It can generate several output files all with a common prefix and different endings:

- per-base depth (`{prefix}.per-base.bed.gz`),
- mean per-window depth given a window size (`{prefix}.regions.bed.gz`, if a BED file provided with `--by`),
- mean per-region given a BED file of regions (`{prefix}.regions.bed.gz`, if a window size provided with `--by`),
- a distribution of proportion of bases covered at or above a given threshhold for each chromosome and genome-wide (`{prefix}.mosdepth.global.dist.txt` and `{prefix}.mosdepth.region.dist.txt`),
- quantized output that merges adjacent bases as long as they fall in the same coverage bins (`{prefix}.quantized.bed.gz`),
- threshold output to indicate how many bases in each region are covered at the given thresholds (`{prefix}.thresholds.bed.gz`)

The MultiQC module plots coverage distributions from 2 kinds of outputs:

- `{prefix}.mosdepth.region.dist.txt`
- `{prefix}.mosdepth.global.dist.txt`

Plotting 2 plots for each:

- Distribution of the number of locations in the genome with a given depth of coverage.
- Average coverage per contig/chromosome.
