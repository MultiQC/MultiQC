---
Name: InterOp
URL: http://illumina.github.io/interop/index.html
Description: >
    The Illumina InterOp libraries are a set of common routines used for reading and writing InterOp metric files. These metric files are binary files produced during a run providing detailed statistics about a run. In a few cases, the metric files are produced after a run during secondary analysis (index metrics) or for faster display of a subset of the original data (collapsed quality scores).
---

This module parses the output from the InterOp Summary executable and creates a table view. The aim is to replicate the `Run & Lane Metrics` table from the [Illumina Basespace](https://basespace.illumina.com) interface. The executable used can easily be installed from the BioConda channel using `conda install -c bioconda illumina-interop`.