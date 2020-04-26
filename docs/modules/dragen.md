---
Name: DRAGEN
URL: https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html
Description: >
    Illumina Bio-IT Platform that uses FPGA for secondary NGS analysis.
---

[Illumina DRAGEN](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html)
is a Bio-IT Platform that provides ultra-rapid secondary analysis of sequencing data using field-programmable 
gate array technology (FPGA).

DRAGEN has a number of differrent pipelines and outputs, including base calling, DNA and RNA alignment,
post-alignment processing and variant calling, covering virtually all stages of typical NGS data processing.
For each stage, it generates QC files with metrics resembling those of samtools-stats, mosdepth, bcftools-stats 
and alike. This MultiQC module supports some of the output but not all.

- `<output prefix>.wgs_fine_hist_<tumor|normal>.csv        - coverage distribution and cumulative coverage plots`
- `<output prefix>.mapping_metrics.csv                     - general stats table, a dedicated table, and a few barplots`
- `<output prefix>.wgs_coverage_metrics_<tumor|normal>.csv - general stats table and a dedicated table`
- `<output prefix>.wgs_contig_mean_cov_<tumor|normal>.csv  - a histogram like in mosdepth, with each chrom as a category on X axis, plus a category for autosomal chromosomes average`
- `<output prefix>.fragment_length_hist.csv                - a histogram plot`
- `<output prefix>.ploidy_estimation_metrics.csv           - add just Ploidy estimation into the general stats table`
- `<output prefix>.vc_metrics.csv                          - a dedicated table and the total number of Variants into the general stats table`
- `<output prefix>.time_metrics.csv                        - a beeswarm plot with dots for time each stage for each sample has finished`

Each QC output adds a section into the report if a corresponding QC file is found.
