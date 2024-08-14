---
name: Sickle
urls: ["https://github.com/najoshi/sickle"]
summary: >
  A windowed adaptive trimming tool for FASTQ files using quality
---

The `stdout` can be captured by directing it to a file e.g. `sickle command 2> sickle_out.log`

The module generates the sample names based on the filenames.

### File search patterns

```yaml
sickle:
  contents_re: 'FastQ \w*\s?records kept: .*'
  num_lines: 2
```
