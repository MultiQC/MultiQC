---
name: Cluster Flow
urls: ["http://clusterflow.io"]
summary: >
  Simple and flexible bioinformatics pipeline tool
---

The module for Cluster Flow parses `*_clusterflow.txt` logs
and finds consensus commands executed by modules in each pipeline run.

The Cluster Flow `*.run` files are also parsed and pipeline information
shown (some basic statistics plus the pipeline steps / params used).

### File search patterns

```yaml
clusterflow/logs:
  fn: "*_clusterFlow.txt"
  shared: true
clusterflow/runfiles:
  contents: Cluster Flow Run File
  fn: "*.run"
  num_lines: 2
```