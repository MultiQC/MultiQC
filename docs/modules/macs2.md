---
Name: MACS2
URL: https://github.com/taoliu/MACS
Description: >
    MACS2 identifies transcription factor binding sites in ChIP-seq data.
---

MACS2 _(Model-based Analysis of ChIP-Seq)_ is a tool for identifying transcript
factor binding sites. MACS captures the influence of genome complexity to
evaluate the significance of enriched ChIP regions.

The MACS2 MultiQC module reads the header of the `*_peaks.xls` results files
and prints the redundancy rates in the General Statistics table. Numerous
additional values are parsed and saved to `multiqc_data/multiqc_macs2.txt`.
