---
name: ngs-bits
urls: ["https://github.com/imgag/ngs-bits"]
summary: >
  Calculating statistics from FASTQ, BAM, and VCF
---

The ngs-bits module parses XML output generated for several tools in the ngs-bits collection:

- [ReadQC](https://github.com/imgag/ngs-bits/blob/master/doc/tools/ReadQC.md) for statistics on FASTQ files,
- [MappingQC](https://github.com/imgag/ngs-bits/blob/master/doc/tools/MappingQC.md) for statistics on BAM files

### File search patterns

```yaml
ngsbits/mappingqc:
  - contents: MappingQC
    fn: "*.qcML"
    num_lines: 20
ngsbits/readqc:
  - contents: ReadQC
    fn: "*.qcML"
    num_lines: 20
  - contents: SeqPurge
    fn: "*.qcML"
    num_lines: 20
```
