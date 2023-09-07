---
title: Introduction
description: How to install MultiQC on your system
---

# Introduction

MultiQC is a reporting tool that parses results and statistics from bioinformatics tool outputs, such as log files and console outputs.
It helps to summarise experiments containing multiple samples and multiple analysis steps.
It's designed to be placed at the end of pipelines or to be run manually when you've finished running your tools.

![MultiQC Overview](../../images/multiqc_overview.excalidraw.svg)

:::note
MultiQC doesn't _do_ any analysis for you - it just finds results from other tools that you have already run and generates nice reports.
:::

When you launch MultiQC, it recursively searches through any provided file paths and finds files that it recognises. It parses relevant information from these and generates a single stand-alone HTML report file.

In addition to the HTML report, MultiQC generates a directory of parsed data files with consistent data structure. This can be useful for further downstream analysis.
