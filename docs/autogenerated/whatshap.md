---
name: WhatsHap
urls: ["https://whatshap.readthedocs.io/"]
summary: >
  Phasing genomic variants using DNA reads (aka read-based phasing, or haplotype assembly)
---

The module is currently restricted to the output from `whatshap stats --tsv`.

### File search patterns

```yaml
whatshap/stats:
  contents: "#sample\tchromosome\tfile_name\tvariants\tphased\tunphased\tsingletons"
  num_lines: 1
```
