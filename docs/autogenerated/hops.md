---
name: HOPS
urls: ['https://github.com/rhuebler/HOPS/']
summary: >
  Ancient DNA characteristics screening tool of output from the metagenomic aligner MALT
---

This module takes the JSON output of the HOPS postprocessing R script (version >= 0.34) to recreate the
possible positives heatmap, with the heat intensity representing the number of 'ancient DNA characteristics'
categories (small edit distance, damage, both edit distance and aDNA damage) that a particular taxon has.

### File search patterns

```yaml
hops:
  fn: heatmap_overview_Wevid.json
```