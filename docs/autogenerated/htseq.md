---
name: HTSeq Count
urls: ["https://htseq.readthedocs.io/en/master/htseqcount.html"]
summary: >
  Part of the HTSeq package: counts reads covering specified genomic features
---

HTSeq is a general purpose Python package that provides infrastructure to
process data from high-throughput sequencing assays. `htseq-count` is a tool
that is part of the main HTSeq package - it takes a file with aligned sequencing
reads, plus a list of genomic features and counts how many reads map to each feature.

### File search patterns

```yaml
htseq:
  - contents_re: ^feature\tcount$
    num_lines: 1
    shared: true
  - contents_re: ^\w+.*\t\d+$
    num_lines: 1
    shared: true
```