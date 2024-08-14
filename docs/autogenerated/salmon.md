---
name: Salmon
urls: ['https://combine-lab.github.io/salmon/']
summary: >
  Quantifies expression of transcripts using RNA-seq data
---

The Salmon module parses `meta_info.json`, `lib_format_counts.json` and `flenDist.txt` files, if found.

:::note
Note that `meta_info.json` must be within a directory called either `aux_info` or `aux` and will be ignored
otherwise.
:::

### File search patterns

```yaml
salmon/fld:
  fn: flenDist.txt
salmon/lfc:
  fn: lib_format_counts.json
salmon/meta:
  contents: salmon_version
  fn: meta_info.json
  max_filesize: 50000
  num_lines: 10
```