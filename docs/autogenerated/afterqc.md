---
name: AfterQC
urls: ["https://github.com/OpenGene/AfterQC"]
summary: >
  Automatic filtering, trimming, error removing, and quality control for FastQ data
---

AfterQC goes through all FastQ files in a folder and outputs three folders: good, bad and QC folders,
which contains good reads, bad reads and the QC results of each fastq file/pair.

### File search patterns

```yaml
afterqc:
  contents: allow_mismatch_in_poly
  fn: "*.json"
  num_lines: 10000
```
