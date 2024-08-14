---
name: Conpair
urls: ["https://github.com/nygenome/Conpair"]
summary: >
  Estimates concordance and contamination for tumor–normal pairs
---

Useful for tumor-normal studies. Performs concordance verification (= samples coming from the same individual), and cross-individual contamination level estimation in WGS and WES sequencing experiments

### File search patterns

```yaml
conpair/concordance:
  contents: markers (coverage per marker threshold
  num_lines: 3
conpair/contamination:
  contents: "Tumor sample contamination level: "
  num_lines: 3
```
