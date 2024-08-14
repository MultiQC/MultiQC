---
name: HUMID
urls: ["https://github.com/jfjlaros/HUMID"]
summary: >
  Reference-free tool to quickly remove duplicates from FastQ files, with or without UMIs
---

### File search patterns

```yaml
humid/clusters:
  contents_re: "[0-9]+ [0-9]+"
  fn: clusters.dat
  num_lines: 1
humid/counts:
  contents_re: "[0-9]+ [0-9]+"
  fn: counts.dat
  num_lines: 1
humid/neighbours:
  contents_re: "[0-9]+ [0-9]+"
  fn: neigh.dat
  num_lines: 1
humid/stats:
  contents: "total: "
  fn: stats.dat
  num_lines: 1
```
