---
title: Preseq
displayed_sidebar: multiqcSidebar
description: >
  Estimates library complexity, showing how many additional unique reads are sequenced for increasing total read count.
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/preseq/preseq.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

:::note
Estimates library complexity, showing how many additional unique reads are sequenced for increasing total read count.

[http://smithlabresearch.org/software/preseq/](http://smithlabresearch.org/software/preseq/)
:::

A shallow curve indicates complexity saturation. The dashed line shows a perfectly complex library where total reads = unique reads.

When `preseq lc_extrap` is run with the default parameters, the extrapolation points
reach 10 billion molecules making the plot difficult to interpret in most scenarios.
It also includes a lot of data in the reports, which can unnecessarily inflate report
file sizes. To avoid this, MultiQC trims back the x-axis until each dataset
shows 80% of its maximum y-value (unique molecules).

To disable this feature and show all the data, add the following to your
[MultiQC configuration](https://docs.seqera.io/multiqc/getting_started/config):

```yaml
preseq:
  notrim: true
```

#### Using coverage instead of read counts

Preseq reports its numbers as "Molecule counts". This isn't always very intuitive,
and it's often easier to talk about sequencing depth in terms of coverage.
You can plot the estimated coverage instead by specifying the reference genome or target size,
and the read length in your [MultiQC configuration](https://docs.seqera.io/multiqc/getting_started/config):

```yaml
preseq:
  genome_size: 3049315783
  read_length: 300
```

These parameters make the script take every molecule count and divide it by
(genome_size / read_length).

MultiQC comes with effective genome size presets for Human and Mouse, so you can
provide the genome build name instead, like this: `genome_size: hg38_genome`. The
following values are supported: `hg19_genome`, `hg38_genome`, `mm10_genome`.

When the genome and read sizes are provided, MultiQC will plot the molecule counts
on the X axis ("total" data) and coverages on the Y axis ("unique" data).
However, you can customize what to plot on each axis (counts or coverage), e.g.:

```yaml
preseq:
  x_axis: counts
  y_axis: coverage
```

#### Plotting externally calculated read counts

To mark on the plot the read counts calculated externally from BAM or fastq files,
create a file with `preseq_real_counts` in the filename and place it with your analysis files.
It should be space or tab delimited with 2 or 3 columns (column 1 = preseq file name,
column 2 = real read count, optional column 3 = real unique read count). For example:

```
Sample_1.preseq.txt 3638261 3638011
Sample_2.preseq.txt 1592394 1592133
[...]
```

You can generate a line for such a file using samtools:

```bash
echo "Sample_1.preseq.txt "$(samtools view -c -F 4 Sample_1.bam)" "$(samtools view -c -F 1028 Sample_1.bam)
```

### File search patterns

```yaml
preseq:
  - contents: EXPECTED_DISTINCT
    num_lines: 2
  - contents: distinct_reads
    num_lines: 2
preseq/real_counts:
  fn: "*preseq_real_counts*"
```
