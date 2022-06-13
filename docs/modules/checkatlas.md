---
Name: CheckAtlas
URL: https://github.com/becavin-lab/checkatlas
Description: >
A one-liner tool for quality control of your single-cell atlases.
---

CheckAtlas is a one liner tool to check the quality of your single-cell atlases. For every atlas, it produces the
quality control tables and figures which can be then processed by multiqc. CheckAtlas is able to load Scanpy, Seurat,
and CellRanger files.

Multiqc will parse all this different tables and figures:

- `myatlas_checkatlas_summ.tsv` - Summary tables with general information on atlases
- `myatlas_checkatlas_adata.tsv`- Table with all scanpy adata features
- `mydata_checkatlas_qc.png` - Quality control figures for every atlas
- `mydata_checkatlas_umap.png` - All UMAP reductions found in the atlas
- `mydata_checkatlas_tsne.png` - All t-SNE reductions found in the atlas
- `mydata_checkatlas_mcluster.tsv` - Table with cluster metrics calculated for every atlas
- `mydata_checkatlas_mannot.tsv` - Table with annotation metrics calculated for every atlas
- `mydata_checkatlas_mdimred.tsv` - Table with dimensionality reduction metrics calculated for every atlas
- `mydata_checkatlas_mspecificity.tsv` - Table with specificity metrics calculated for every atlas

`multiqc --cl-config "ignore_images: false"` option should be added to get all checkatlas output in multiqc.
