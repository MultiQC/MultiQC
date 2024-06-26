---
name: Kraken
url: https://ccb.jhu.edu/software/kraken2/
description: >
  is a taxonomic classification tool that uses exact k-mer matches to find
  the lowest common ancestor (LCA) of a given sequence.
---

The MultiQC module supports outputs from both Kraken and Kraken 2.

It works with report files generated using the `--report` flag, that look like the following:

```ts
11.66	98148	98148	U	0	unclassified
88.34	743870	996	-	1	root
88.22	742867	0	-	131567	  cellular organisms
88.22	742866	2071	D	2	    Bacteria
87.95	740514	2914	P	1239	      Firmicutes
```

A bar graph is generated that shows the number of fragments for each sample that
fall into the top-5 categories for each taxa rank. The top categories are calculated
by summing the library percentages across all samples.

The number of top categories to plot can be customized in the config file:

```yaml
kraken:
  top_n: 5
```
