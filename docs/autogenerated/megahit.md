---
name: MEGAHIT
urls: ["https://github.com/voutcn/megahit"]
summary: >
  NGS read assembler
---

MultiQC will parse stdout/stderr logs from MEGAHIT runs. The sample name is taken from the file
name (e.g. `sample1.log` will yield a sample name of `sample1`).

### File search patterns

```yaml
megahit:
  contents: " - MEGAHIT v"
  num_lines: 5
```
