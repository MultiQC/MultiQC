---
Name: mosdepth
URL: https://github.com/brentp/mosdepth
Description: >
    Mosdepth performs fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing.
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

Using "region" if exists, otherwise "global". Plotting 3 figures:

- Distribution of the number of locations in the genome with a given depth of coverage.
- Absoulute number of locations in the genome with a given depth of coverage.
- Average coverage per contig/chromosome.

Also plotting the percentage of the genome covered at a threshold in the General Stats section. 
The default thresholds are 1, 5, 10, 30, 50, which can be customised in the config as follows:

```yaml
mosdepth_config:
    general_stats_coverage:
        - 10
        - 20
        - 40
        - 200
        - 30000
```

You can also specify which columns would be hidden when the report loads (by default, all values are hidden except 30X):

```yaml
    general_stats_coverage_hidden:
        - 10
        - 20
        - 200
```
