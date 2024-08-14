---
name: Bakta
urls: ["https://github.com/oschwengers/bakta"]
summary: >
  Rapid & standardized annotation of bacterial genomes, MAGs & plasmids
---

The module analyses summary results from the Bakta annotation pipeline for bacterial genomes. The
summary text file used is included in the Bakta output since v1.3.0. The MultiQC module was written for
the output of v1.7.0.

### File search patterns

```yaml
bakta:
  contents: "Bakta:"
  fn: "*.txt"
```
