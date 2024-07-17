---
name: GoPeaks
url: https://github.com/maxsonBraunLab/gopeaks
description: Calls peaks in CUT&TAG/CUT&RUN datasets
---

Gopeaks uses a binomial distribution to model the read counts in sliding windows across
the genome and calculate peak regions that are enriched over the background.

The module recognizes files with the `*_gopeaks.json` suffix (which is the default behavior), and will report
the number of peaks called per sample via the general table and the bar plot.
