---
Name: Sentieon
URL: https://www.sentieon.com/products/
Description: >
    Sentieon-dnaseq produces many outputs. This module deals with 3 Picard
    equivalents which do not transfer well to MultiQC. The code for each script
    is split into its own file and adds a section to the module output if
    logs are found.
---

The Sentieon module parses output from the Sentieon dna-seq suite of tools,
which themselves are implementations of certain Picard metrics
[Picard](http://broadinstitute.github.io/picard/),
a set of Java command line tools for manipulating high-throughput
sequencing data.

Supported commands:

* `InsertSizeMetrics`
* `GcBiasMetrics`
* `AlignmentSummaryMetrics`

