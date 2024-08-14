---
name: RSEM
urls: ["https://deweylab.github.io/RSEM/"]
summary: >
  Estimates gene and isoform expression levels from RNA-Seq data
---

Supported scripts:

- `rsem-calculate-expression`

This module search for the file `.cnt` created by RSEM into directory named `PREFIX.stat`

### File search patterns

```yaml
rsem:
  fn: "*.cnt"
```
