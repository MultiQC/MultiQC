---
Name: GffCompare
URL: https://ccb.jhu.edu/software/stringtie/gffcompare.shtml
Description: A tool to compare, merge and annotate one or more GFF files with a reference annotation in GFF format.
---

The program `gffcompare` can be used to compare, merge, annotate and estimate accuracy
of one or more GFF files (the "query" files), when compared with a reference annotation (also provided as GFF).

The _Sensitivity / Precision_ values are displayed in a single plot,
different loci levels can be switched by choosing a different dataset.

> **NB:** Please use `gffcompare` only with single samples.
> Multi-Sample comparisons are not correctly rendered by this MultiQC module.

Note that exported data in `multiqc_data/multiqc_gffcompare.{tsv,yaml,json}` only works when
exporting with YAML or JSON - the default `.tsv` output will not contain any data.
Please use `-k yaml` or `-k json` to export in a structured format.
It is hoped to refactor this code in a future release - please submit a PR if you are interested.
