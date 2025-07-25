---
title: Bustools
displayed_sidebar: multiqcSidebar
description: >
  Tools for BUS files - a file format for single-cell RNA-seq data designed to facilitate the development of modular workflows for data processing.
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/bustools/bustools.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

:::note
Tools for BUS files - a file format for single-cell RNA-seq data designed to facilitate the development of modular workflows for data processing.

[https://bustools.github.io/](https://bustools.github.io/)
:::

This module parses the report generated by Bustools `inspect`, and expects the file to be
named `inspect.json`. This is the default naming pattern when you make use of the
[kallisto-bustools wrapper](https://www.kallistobus.tools/). The sample name is set as the name
of the directory containing the file. For the rest, this module should work in exactly the same
way as all other MultiQC modules.

### File search patterns

```yaml
bustools:
  fn: "*inspect.json"
```
