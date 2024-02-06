---
name: metaDMG
url: https://github.com/metaDMG-dev/metaDMG-cpp
description: >
  metaDMG-cpp is a fast and efficient method for assigning reads to LCA and to estimate damage rates in ancient DNA data.
---

It relies on a commonly used alignment files formats bam/sam/sam.gz and can calculate
degree of damage from read data mapped against a single as well as multiple genomes.
It is especially relevant for data where data have been mapped against multiple
reference genomes or to speed up analysis for damage estimation for a single genome.

This module can read both plain text and gzipped files. If the latter, to allow reading
the gz archives, run with `ignore_images: false` in the config, e.g.:

```
multiqc . --cl-config 'ignore_images: false'
```
