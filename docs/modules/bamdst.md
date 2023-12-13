---
name: Bamdst
url: https://https://github.com/shiquan/bamdst
description: >
  Bamdst is a lightweight tool to stat the depth coverage of target regions of bam file(s).
---

The MultiQC module reads data from two types of Bamdst logs:

- `coverage.report`: used to build a table with coverage statistics. The sample name is read from this file.
- `chromosomes.report`: if this file is found in the same directory as the file above, additionally a per-contig coverage plot will be generated.

For the per-contig coverage plot, you can include and exclude contigs based on name or pattern.

For example, you could add the following to your MultiQC config file:

```yaml
bamdst:
  include_contigs:
    - "chr*"
  exclude_contigs:
    - "*_alt"
    - "*_decoy"
    - "*_random"
    - "chrUn*"
    - "HLA*"
    - "chrM"
    - "chrEBV"
```

Note that exclusion supersedes inclusion for the contig filters.

To additionally avoid cluttering the plot, MultiQC can exclude contigs with a low relative coverage.

```yaml
bamdst:
  # Should be a fraction, e.g. 0.001 (exclude contigs with 0.1% coverage of sum of
  # coverages across all contigs)
  perchrom_fraction_cutoff: 0.001
```

If you want to see what is being excluded, you can set `show_excluded_debug_logs` to `True`:

```yaml
bamdst:
  show_excluded_debug_logs: True
```

This will then print a debug log message (use `multiqc -v`) for each excluded contig.
This is disabled by default as there can be very many in some cases.
