---
name: Motus
urls: ["https://motu-tool.org/"]
summary: >
  Microbial profiling through marker gene (MG)-based operational taxonomic units (mOTUs)
---

The module takes as input in the stdout of `mOTUs profile`, and provides summary statistics on various steps of the pipeline.

### File search patterns

```yaml
motus:
  contents: Reads are aligned (by BWA) to marker gene sequences in the reference database
  num_lines: 2
```
