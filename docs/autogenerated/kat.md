---
name: K-mer Analysis Toolkit
urls: ["https://github.com/TGAC/KAT"]
summary: >
  Analyses sequencing data via its k-mer spectra
---

The KAT multiqc module interprets output from KAT distribution analysis json files, which typically
contain information such as estimated genome size and heterozygosity rates from your k-mer spectra.

### File search patterns

```yaml
kat:
  fn: "*.dist_analysis.json"
```
