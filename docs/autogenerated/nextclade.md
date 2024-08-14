---
name: Nextclade
urls: ["https://github.com/nextstrain/nextclade"]
summary: >
  Viral genome alignment, clade assignment, mutation calling, and quality checks
---

Nextclade assigns input sequences to SARS-Cov-2 clades based on differences between the input sequences
and [Nextstrain](https://nextstrain.org/) reference sequences. In addition, it judges the validity of
the samples by performing several quality control checks on the input sequences.

### File search patterns

```yaml
nextclade:
  contents: seqName;clade;
  num_lines: 1
```
