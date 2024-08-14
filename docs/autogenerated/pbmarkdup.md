---
name: pbmarkdup
urls: ['https://github.com/PacificBiosciences/pbmarkdup']
summary: >
  Takes one or multiple sequencing chips of an amplified libray as HiFi reads and marks or removes duplicates
---

The module adds the **% Unique Molecules** and **%Duplicate Reads** (hidden) to the General Statistics
table.

### File search patterns

```yaml
pbmarkdup:
  contents_re: LIBRARY +READS +UNIQUE MOLECULES +DUPLICATE READS
  num_lines: 5
```