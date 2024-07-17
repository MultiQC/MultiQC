---
name: DRAGEN
url: https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html
description: >
  Illumina Bio-IT Platform that uses FPGA for secondary analysis of sequencing data
---

DRAGEN has a number of different pipelines and outputs, including base calling, DNA and RNA alignment,
post-alignment processing and variant calling, covering virtually all stages of typical NGS data processing.
For each stage, it generates QC files with metrics resembling those of samtools-stats, mosdepth, bcftools-stats
and alike. This MultiQC module supports some of the output but not all. Contributions are welcome!

- `<output prefix>.wgs_fine_hist_<tumor|normal>.csv`
  - Coverage distribution and cumulative coverage plots
- `<output prefix>.mapping_metrics.csv`
  - General stats table, a dedicated table, and a few barplots
- `<output prefix>.wgs_coverage_metrics_<tumor|normal>.csv`
  - General stats table and a dedicated table
- `<output prefix>.qc-coverage-region-<1|2|3>_coverage_metrics.csv`
  - General stats table and a dedicated table
- `<output prefix>.wgs_contig_mean_cov_<tumor|normal>.csv`
  - A histogram like in mosdepth, with each chrom as a category on X axis, plus a category for autosomal chromosomes average
- `<output prefix>.fragment_length_hist.csv`
  - A histogram plot
- `<output prefix>.ploidy_estimation_metrics.csv`
  - Add just Ploidy estimation into the general stats table
- `<output prefix>.vc_metrics.csv`
  - A dedicated table and the total number of Variants into the general stats table
- `<output prefix>.gc_metrics.csv`
  - A histogram and summary statistics table on GC content metrics
- `<output prefix>.trimmer_metrics.csv`
  - A summary table of tirmmer metrics
- `<output prefix>.time_metrics.metrics`
  - A bar graph of the total run time and a breakdown of the run time of each individual step
- `<output prefix>.quant.metrics.csv`
  - A bar graph of RNA fragments
- `<output prefix>.quant.transcript_coverage.txt`
  - A line plot of average coverage along RNA transcripts
- `<output prefix>.scRNA.metrics.csv` or `<output prefix>.scRNA_metrics.csv`
  - Summary table for single-cell RNA metrics
- `<output prefix>.scATAC.metrics.csv` or `<output prefix>.scATAC_metrics.csv`
  - Summary table for single-cell ATAC metrics

The code is structured in a way so every mix-in parses one type of QC file that DRAGEN generates
(e.g. _.mapping_metrics.csv, _.wgs_fine_hist_normal.csv, etc.). If a corresponding file is found, a mix-in adds
a section into the report.

DRAGEN can be treated as a fast aligner with additional features on top, as users will unlikely use any
features without enabling DRAGEN mapping. So we will treat this module as an alignment tool module and
place it accordingly in the module_order list, in docs, etc.
