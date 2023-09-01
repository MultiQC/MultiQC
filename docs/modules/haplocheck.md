---
name: Haplocheck
url: https://github.com/genepi/haplocheck
description: >
  Detects in-sample contamination in mtDNA or WGS sequencing studies by analyzing the mitchondrial content.
---

By parsing the `[sample_name].raw.txt` files, the Haplocheck module generates a table and plots in the MultiQC report.

Please note that it requires the raw output. Use the `--raw` parameter when calling the program.

`haplocheck --raw --out $sample_id $vcf`

MultiQC has been tested with output from Haplocheck 1.3.3.
