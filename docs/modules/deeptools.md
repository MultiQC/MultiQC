---
Name: deepTools
URL: http://deeptools.readthedocs.io
Description: >
    Tools to process and analyze deep sequencing data.
---

deepTools addresses the challenge of handling the large amounts of data that are now routinely generated from DNA sequencing centers. deepTools contains useful modules to process the mapped reads data for multiple quality checks, creating **normalized coverage files** in standard bedGraph and bigWig file formats, that allow comparison between different files (for example, treatment and control). Finally, using such normalized and standardized files, deepTools can create many publication-ready **visualizations** to identify enrichments and for functional annotations of the genome.

The MultiQC module for deepTools parses a number of the text files that deepTools can produce. In particular, the following are supported:

 - `bamPEFragmentSize --table`
 - `estimateReadFiltering`
 - `plotCoverage ---outRawCounts` (as well as the content written normally to the console)
 - `plotEnrichment --outRawCounts`
 - `plotFingerprint --outQualityMetrics --outRawCounts`

Please be aware that some tools (namely, `plotFingerprint --outRawCounts` and `plotCoverage --outRawCounts`) are only supported as of deepTools version 2.6. For earlier output from `plotCoverage --outRawCounts`, you can use `#'chr'	'start'	'end'	` in `utils/search_patterns.yaml` (see [here](http://multiqc.info/docs/#module-search-patterns) for more details). Also for these types of files, you may need to increase the maximum file size supported by MultiQC (`log_filesize_limit` in the MultiQC configuration file). You can find details regarding the configuration file location [here](http://multiqc.info/docs/#configuring-multiqc).

Note that sample names are parsed from the text files themselves, they are not derived from file names.
