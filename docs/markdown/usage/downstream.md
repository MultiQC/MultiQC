---
title: Downstream analysis
description: How to use MultiQC raw data outputs
---

# Downstream analysis

Whilst MultiQC is typically used as a final reporting step in an analysis, it can also be used as an intermediate in your analysis.

MultiQC saves a directory of machine-readable outputs called `multiqc_data/`. In here there are files from each module and table, as well as a verbose `multiqc.log` file and a `multiqc_data.json` file that contains just about everything.

Most of these files are tab-separated `.tsv` files by default, but you can choose to have them as JSON, YAML if you prefer with the `-k`/`--data-format` flag or the `data_format` option in a config file.

These files can be useful as MultiQC essentially standardises the outputs from a lot of different tools.
Typical usage of MultiQC outputs could be filtering of large datasets (eg. single-cell analysis) or trend-monitoring of repeated runs.

Below are a few tools that are specifically designed to work with MultiQC.
They are not created by or endorsed by the MultiQC author but may be helpful for your research.

## TidyMultiqc

- Homepage: [https://CRAN.R-project.org/package=TidyMultiqc](https://CRAN.R-project.org/package=TidyMultiqc)
- Source: [https://github.com/TMiguelT/TidyMultiqc](https://github.com/TMiguelT/TidyMultiqc)

Provides the means to convert `multiqc_data.json` files into `tidy` data frames for downstream analysis in R.

This analysis might involve cohort analysis, quality control visualisation, change-point detection, statistical process control, clustering, or any other type of quality analysis.

## MegaQC

- Homepage: [https://megaqc.info](https://megaqc.info)
- Source: [https://github.com/ewels/MegaQC](https://github.com/ewels/MegaQC)

Started off by MultiQC author [@ewels](https://github.com/ewels/) this project has had further development by a team of several contributors. It is functional but still has several parts of its codebase that have never quite been finished.

MegaQC imports data from multiple MultiQC runs and provides an interface to explore this with an interactive web server using a database backend.
It can plot data over time, across runs and even has an interactive dashboard builder.
It's useful for anyone who wants to monitor MultiQC statistics (eg. clinical labs) or work interactively with large datasets (eg. single cell analysis).

## ChronQC

- Docs: [https://chronqc.readthedocs.io](https://chronqc.readthedocs.io)
- Source: [https://github.com/nilesh-tawari/ChronQC](https://github.com/nilesh-tawari/ChronQC)

ChronQC is a quality control (QC) tracking system for clinical implementation of next-generation sequencing (NGS). ChronQC generates time series plots for various QC metrics, which allows comparison of the current run to historical runs. ChronQC has multiple features for tracking QC data including Westgard rules for clinical validity, laboratory-defined thresholds, and historical observations within a specified period. Users can record their notes and corrective actions directly onto the plots for long-term recordkeeping.
