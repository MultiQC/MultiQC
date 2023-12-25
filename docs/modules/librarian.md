---
name: Librarian
url: https://github.com/DesmondWillowbrook/Librarian
description: >
  A tool to predict the sequencing library type from the base
  composition of a supplied FastQ file.
---

Librarian reads from high throughput sequencing experiments show base compositions that are characteristic for their library type.
For example, data from RNA-seq and WGBS-seq libraries show markedly different distributions of G, A, C and T across the reads.

Librarian makes use of different composition signatures for library quality control: Test library compositions are extracted and compared against previously published data sets from mouse and human.

This MultiQC module generates the _Prediction Plot_ showing the likelihood that samples are a given library type.

### General Stats

The module can also show the most likely library type in the General Statistics table, however this is disabled by default.
This is because several library types are very similar to each other and can come out as a mix.
It's often misleading to show only the top one (even if it has a low score), but very clear when looking at the full heatmap.

If you really want to show the most likely library type, you can enable this in your MultiQC config file:

```yaml
librarian:
  show_general_stats: true
```
