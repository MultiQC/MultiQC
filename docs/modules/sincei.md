---
name: sincei
url: http://sincei.readthedocs.io
description: >
  A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data
---

sincei (short for **Single Cell Informatics**) is a command-line toolkit for exploration of single-cell epigenomics data. It accommodates data from a wide range of single-cell protocols, such as droplet-based (10x genomics) and plate-based protocols, gene expression (scRNA-seq) and epigenomics (scATAC, scCUTnTAG, scBS-seq) protocols. sincei can be used for quality control of these datasets directly from BAM files (read-level QC),  as well as after signal aggregation (count-level). Additional functionalities include filtering, dimentionality reduction, and plotting of single-cell data.

The MultiQC module for sincei parses a number of the text files produced by sincei quality control modules. In particular, files produced by:

- `scFilterStats` (default output file)
- `scCountQC --outMetrics` (currently the `cell-level` metrics are supported)


**Note**
  - Each cell (marked by `Cell_ID` in sincei) is considered a "sample" by multiQC.
  - sample/cell names are parsed from the text files themselves, they are not derived from file names.
