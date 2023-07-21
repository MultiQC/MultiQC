---
name: CheckAtlas
url: https://github.com/becavin-lab/checkatlas
description: >
  A one-liner tool for quality control of your single-cell atlases.
---

CheckAtlas is a one liner tool to check the quality of your single-cell atlases. For every atlas, it produces the
quality control tables and figures which can be then processed by multiqc. CheckAtlas is able to load Scanpy, Seurat,
and CellRanger files.

Multiqc will parse all these different tables:

- `summary/sample_name.tsv` - Summary tables with general information on atlases
- `adata/sample_name.tsv`- Table with all scanpy adata features
- `qc/sample_name.tsv` - Quality control tables for every atlas
- `cluster/sample_name.tsv` - Table with cluster metrics calculated for every atlas
- `annot/sample_name.tsv` - Table with annotation metrics calculated for every atlas
- `dimred/sample_name.tsv` - Table with dimensionality reduction metrics calculated for every atlas
