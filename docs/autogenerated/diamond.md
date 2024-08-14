---
name: DIAMOND
urls: ["https://github.com/bbuchfink/diamond"]
summary: >
  Sequence aligner for protein and translated DNA searches, a drop-in replacement for the NCBI BLAST
---

Key features are:

- Pairwise alignment of proteins and translated DNA at 100x-10,000x speed of BLAST.
- Frameshift alignments for long read analysis.
- Low resource requirements and suitable for running on standard desktops or laptops.
- Various output formats, including BLAST pairwise, tabular and XML, as well as taxonomic classification.

The module takes summary statistics from the `diamond.log` file (`--log` option). It parses and reports
the number of sequences aligned and displays them in the General Stats table.

### File search patterns

```yaml
diamond:
  fn: diamond.log
```
