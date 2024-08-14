---
name: MinIONQC
urls: ["https://github.com/roblanf/minion_qc"]
summary: >
  Quality control for ONT (Oxford Nanopore) long reads
extra_description: >
  It uses the `sequencing_summary.txt` files produced by ONT (Oxford Nanopore Technologies)
  long-read base-callers to perform QC on the reads. It allows quick-and-easy comparison of data from
  multiple flowcells
---

The module parses data in the `summary.yaml` MinIONQC output files.

### File search patterns

```yaml
minionqc:
  contents: total.gigabases
  fn: summary.yaml
```
