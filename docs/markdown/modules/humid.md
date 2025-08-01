---
title: HUMID
displayed_sidebar: multiqcSidebar
description: >
  Reference-free tool to quickly remove duplicates from FastQ files, with or without UMIs.
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/humid/humid.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

:::note
Reference-free tool to quickly remove duplicates from FastQ files, with or without UMIs.

[https://github.com/jfjlaros/HUMID](https://github.com/jfjlaros/HUMID)
:::

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
