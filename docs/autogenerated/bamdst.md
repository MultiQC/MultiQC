---
name: Bamdst
urls: ["https://https://github.com/shiquan/bamdst"]
summary: >
  Lightweight tool to stat the depth coverage of target regions of BAM file(s)
---

The module reads data from two types of Bamdst logs:

- `coverage.report`: used to build a table with coverage statistics. The sample name is read from this file.
- `chromosomes.report`: if this file is found in the same directory as the file above, additionally a per-contig coverage plot will be generated. This file must be named exactly this way, with the `.report` extension.

Note that for the sample names, the module will attempt to use the input BAM name
in the header in the `coverage.report` file:

```
## The file was created by bamdst
## Version : 1.0.9
## Files : ST0217_Lg.bam
...
```

However, if the tool was run in a piped manner, the file name will be just `-` or `/dev/stdin`,
and instead MultiQC will fall back to using the log file name `coverage.report`.
Make sure to run MultiQC with `--dirs` if use have multiple samples run in this way,
otherwise MultiQC will only report the first found sample under the name `coverage`.

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
    - "*_fix"
    - "HLA*"
    - "chrUn*"
    - "chrEBV"
    - "chrM"
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

### File search patterns

```yaml
bamdst/coverage:
  contents: "## The file was created by bamdst"
  num_lines: 5
```
