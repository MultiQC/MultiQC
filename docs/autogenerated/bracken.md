---
name: Bracken
urls: ["https://ccb.jhu.edu/software/bracken/"]
summary: >
  Computes the abundance of species in DNA sequences from a metagenomics sample
---

This module works with Bracken output files that resemble Kraken reports. They look like the following:

```tsv
100.00	1188381	0	R	1	root
100.00	1188380	0	R1	131567	  cellular organisms
100.00	1188328	0	D	2	    Bacteria
99.99	1188317	0	P	1224	      Proteobacteria
99.99	1188308	0	C	1236	        Gammaproteobacteria
99.99	1188285	0	O	91347	          Enterobacterales
99.98	1188147	0	F	543	            Enterobacteriaceae
```

The main assumption to tell Bracken reports from Kraken is that the former don't have
an "unassigned" category at the head of the file, and instead start with "root".

A bar graph is generated identical to that of Kraken.

### File search patterns

```yaml
bracken:
  contents_re: ^(\d{1,3}\.\d{1,2})\t(\d+)\t(\d+)\t((\d+)\t(\d+)\t)?([URDKPCOFGS-]\d{0,2})\t(\d+)(\s+)root
  num_lines: 1
```