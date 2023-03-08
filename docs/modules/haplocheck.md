---
Name: Haplocheck
URL: https://github.com/genepi/haplocheck
Description: >
  Detects in-sample contamination in mtDNA or WGS sequencing studies by analyzing the mitchondrial content.
---

The Haplocheck module parses the `[sample_name].raw.txt` files and creates a table in the MultiQC report.

Please note that it requires the raw output. Use the `--raw` parameter when calling the program.

`haplocheck --raw --out $sample_id $vcf`

MultiQC has been tested with output from Haplocheck 1.3.3.
