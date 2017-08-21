---
Name: THeTA2
URL: http://compbio.cs.brown.edu/projects/theta/
Description: >
    THeTA2 estimates tumour purity and clonal / subclonal copy number.
---

[THeTA2](http://compbio.cs.brown.edu/projects/theta/) (Tumor Heterogeneity Analysis)
is an algorithm that estimates the tumour purity and clonal / subclonal copy number
aberrations directly from high-throughput DNA sequencing data.

The THeTA2 MultiQC module plots the % germline and % tumour subclone for each sample.
Note that each sample can have multiple maximum likelihood solutions - the MultiQC
module plots proportions for the first one in the results file (`*.BEST.results`).
Also note that if there are more than 5 tumour subclones, their percentages are summed.
