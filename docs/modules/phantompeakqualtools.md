---
name: phantompeakqualtools
url: https://www.encodeproject.org/software/phantompeakqualtools/
description: >
  Computes enrichment and quality measures for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data.
---

Used to generate three quality metrics: NSC, RSC, and PBC. The NSC (Normalized strand cross-correlation) and RSC (relative strand cross-correlation) metrics use cross-correlation of stranded read density profiles to measure enrichment independently of peak calling. The PBC (PCR bottleneck coefficient) is an approximate measure of library complexity. PBC is the ratio of (non-redundant, uniquely mappable reads)/(uniquely mappable reads).
