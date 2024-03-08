---
name: Longranger
url: https://www.10xgenomics.com/
description: >
  A set of analysis pipelines that perform sample demultiplexing, barcode processing, alignment, quality control, variant calling, phasing, and structural variant calling.
---

Currently supported Longranger pipelines:

- `wgs`
- `targeted`
- `align`

Usage:

```bash
longranger wgs --fastqs=/path/to/fastq --id=NA12878
multiqc /path/to/NA12878
```

This module will look for the files `_invocation` and `summary.csv` in the the `NA12878` folder, i.e. the output folder of Longranger in this example. The file `summary.csv` is required. If the file `_invocation` is not found the sample will receive a generic name in the MultiQC report (`longranger#1`), instead of `NA12878` or whatever was given by the `--id` parameter.
