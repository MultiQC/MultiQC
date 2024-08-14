---
name: MACS2
urls: ["https://macs3-project.github.io/MACS/"]
summary: >
  Identifies transcription factor binding sites in ChIP-seq data
---

MACS2 _(Model-based Analysis of ChIP-Seq)_ is a tool for identifying transcript
factor binding sites. MACS captures the influence of genome complexity to
evaluate the significance of enriched ChIP regions.

The module reads the `*_peaks.xls` results files and prints the redundancy rates and number of peaks
found in the General Statistics table. Numerous additional values are parsed and saved to
`multiqc_data/multiqc_macs2.txt`.

### File search patterns

```yaml
macs2:
  fn: "*_peaks.xls"
```
