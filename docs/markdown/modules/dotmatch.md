---
title: DotMatch
description: Summarise DotMatch known-target assignment and CRISPR guide-library QC reports.
---

# DotMatch

[DotMatch](https://github.com/dnncha/dotmatch) assigns short DNA sequences to
known targets and keeps ambiguous, unmatched, and invalid outcomes separate.
The MultiQC module summarises the stable report tables emitted by DotMatch
workflow commands.

## Supported files

The module discovers files with the following headers:

- `*sample_qc.tsv` — per-sample read assignment and ambiguity rates;
- `*crispr_qc.summary.tsv` — guide-library coverage, sparsity, Gini, and
  assignment QC.

The module reports a compact set of metrics for visual comparison across
samples. The complete DotMatch TSV and JSON artifacts remain the authoritative
audit record.

## Installation and usage

Install DotMatch from [PyPI](https://pypi.org/project/dotmatch/) or
[Bioconda](https://anaconda.org/bioconda/dotmatch), run the relevant DotMatch
workflow, and then point MultiQC at the output directory:

```bash
multiqc path/to/dotmatch/results
```

The report contains `DotMatch Sample QC` and `DotMatch CRISPR QC` sections when
the corresponding artifacts are present.
