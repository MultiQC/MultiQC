---
title: MEGAHIT
displayed_sidebar: multiqcSidebar
description: >
  NGS read assembler
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/megahit/megahit.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

:::note
NGS read assembler

[https://github.com/voutcn/megahit](https://github.com/voutcn/megahit)
:::

MultiQC will parse stdout/stderr logs from MEGAHIT runs. The sample name is taken from the file
name (e.g. `sample1.log` will yield a sample name of `sample1`).

### File search patterns

```yaml
megahit:
  contents: " - MEGAHIT v"
  num_lines: 5
```
