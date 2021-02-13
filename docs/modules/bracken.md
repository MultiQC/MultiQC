---
Name: Bracken
URL: https://ccb.jhu.edu/software/bracken/
Description: >
  a highly accurate statistical method that computes the abundance of species
  in DNA sequences from a metagenomics sample
---

This module works bracken output files that resemble kraken reports. They look like the following:

```ts
100.00	1188381	0	R	1	root
100.00	1188380	0	R1	131567	  cellular organisms
100.00	1188328	0	D	2	    Bacteria
99.99	1188317	0	P	1224	      Proteobacteria
99.99	1188308	0	C	1236	        Gammaproteobacteria
99.99	1188285	0	O	91347	          Enterobacterales
99.98	1188147	0	F	543	            Enterobacteriaceae
```

A bar graph is generated that shows the number of fragments for each sample that
fall into the top categories for each taxa rank. The top categories are calculated
by summing the library percentages across all samples.
