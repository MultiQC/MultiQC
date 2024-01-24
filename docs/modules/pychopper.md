---
name: Pychopper
url: https://github.com/nanoporetech/pychopper
description: >
  is a tool to identify, orient and trim full-length Nanopore cDNA reads. The tool is also able to rescue fused reads.
---

The MultiQC module parses the pychopper stats file. Pychopper needs to be run with the `-S stats_output` option to create the file. The name of the output file defines the sample name.

The stats file is a three column `tsv` file with the format `category name value`.

Currently only two stats are displayed in MultiQC. Two bargraphs are created for the read classication and the strand orientation of the identified full length transcripts. Additional stats could be included on further request.

The general stats table contains a value that displays the percentage of full length transcripts. This value is calculated from the cumulative length of reads where Pychopper found primers at both ends.
