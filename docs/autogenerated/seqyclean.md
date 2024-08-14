---
name: SeqyClean
urls: ["https://github.com/ibest/seqyclean"]
summary: >
  Filters adapters, vectors, and contaminants while quality trimming
---

SeqyClean is a comprehensive preprocessing software application for NGS reads, that removes noise from FastQ
files to improve de-novo genome assembly and genome mapping.

The module parses the `*SummaryStatistics.tsv` files that results from a SeqyClean cleaning.

### File search patterns

```yaml
seqyclean:
  fn: "*_SummaryStatistics.tsv"
```
