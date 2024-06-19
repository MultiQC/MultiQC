---
name: MetaPhlAn
url: https://github.com/biobakery/MetaPhlAn
description: >
  MetaPhlAn is a computational tool for profiling the composition of microbial communities from metagenomic shotgun sequencing data.
---

The MultiQC module supports outputs from MetaPhlAn, that look like the following:

```tsv
k__Bacteria	2	100.0
k__Bacteria|p__Firmicutes	2|1239	44.30422
k__Bacteria|p__Bacteroidetes	2|976	34.73101
```

A bar graph is generated that shows the relative abundance for each sample that
fall into the top-10 categories for each taxa rank. The top categories are calculated
by summing the relative abundances across all samples.

Any species under the Additional Species column are ignored when making the graphs.

The number of top categories to plot can be customized in the config file:

```yaml
metaphlan:
  top_n: 10
```
