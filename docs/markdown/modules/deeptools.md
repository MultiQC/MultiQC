---
title: deepTools
displayed_sidebar: multiqcSidebar
description: >
  Tools to process and analyze deep sequencing data.
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/deeptools/deeptools.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

:::note
Tools to process and analyze deep sequencing data.

[http://deeptools.readthedocs.io](http://deeptools.readthedocs.io)
:::

deepTools addresses the challenge of handling the large amounts of data that are now routinelygenerated from DNA sequencing centers. deepTools contains useful modules to process the mapped reads data for multiple quality checks, creating **normalized coverage files** in standard bedGraph and bigWig file formats, that allow comparison between different files (for example, treatment and control). Finally, using such normalized and standardized files, deepTools can create many publication-ready **visualizations** to identify enrichments and for functional annotations of the genome.

The module for deepTools parses a number of the text files that deepTools can produce. In particular, the following are supported:

- `bamPEFragmentSize --table`
- `bamPEFragmentSize --outRawFragmentLengths`
- `estimateReadFiltering`
- `plotCoverage ---outRawCounts` (as well as the content written normally to the console)
- `plotEnrichment --outRawCounts`
- `plotFingerprint --outQualityMetrics --outRawCounts`
- `plotPCA --outFileNameData`
- `plotCorrelation --outFileCorMatrix`
- `plotProfile --outFileNameData`

Please be aware that some tools (namely, `plotFingerprint --outRawCounts` and `plotCoverage --outRawCounts`) are only supported as of deepTools version 2.6. For earlier output from `plotCoverage --outRawCounts`, you can use `#'chr' 'start' 'end'` in `search_patterns.yaml` (see [here](https://docs.seqera.io/multiqc/getting_started/config#module-search-patterns) for more details). Also for these types of files, you may need to increase the maximum file size supported by MultiQC (`log_filesize_limit` in the MultiQC configuration file). You can find details regarding the configuration file location [here](https://docs.seqera.io/multiqc/getting_started/config).

Note that sample names are parsed from the text files themselves, they are not derived from file names.

### File search patterns

```yaml
deeptools/bamPEFragmentSizeDistribution:
  contents: "#bamPEFragmentSize"
  num_lines: 1
deeptools/bamPEFragmentSizeTable:
  contents: "\tFrag. Sampled\tFrag. Len. Min.\tFrag. Len. 1st. Qu.\tFrag. Len. Mean\t\
    Frag. Len. Median\tFrag. Len. 3rd Qu."
  num_lines: 1
deeptools/estimateReadFiltering:
  contents: "Sample\tTotal Reads\tMapped Reads\tAlignments in blacklisted regions\t\
    Estimated mapped reads"
  num_lines: 1
deeptools/plotCorrelationData:
  contents: "#plotCorrelation --outFileCorMatrix"
  num_lines: 1
deeptools/plotCoverageOutRawCounts:
  contents: "#plotCoverage --outRawCounts"
  num_lines: 1
deeptools/plotCoverageStdout:
  contents: "sample\tmean\tstd\tmin\t25%\t50%\t75%\tmax"
  num_lines: 1
deeptools/plotEnrichment:
  contents: "file\tfeatureType\tpercent\tfeatureReadCount\ttotalReadCount"
  num_lines: 1
deeptools/plotFingerprintOutQualityMetrics:
  contents: "Sample\tAUC\tSynthetic AUC\tX-intercept\tSynthetic X-intercept\tElbow\
    \ Point\tSynthetic Elbow Point"
  num_lines: 1
deeptools/plotFingerprintOutRawCounts:
  contents: "#plotFingerprint --outRawCounts"
  num_lines: 1
deeptools/plotPCAData:
  contents: "#plotPCA --outFileNameData"
  num_lines: 1
deeptools/plotProfile:
  contents: bin labels
  num_lines: 1
```
