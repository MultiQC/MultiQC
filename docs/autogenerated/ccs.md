---
name: CCS
urls: ["https://github.com/PacificBiosciences/ccs"]
summary: >
  PacBio tool that generates highly accurate single-molecule consensus reads (HiFi Reads)
---

CCS takes multiple subreads of the same SMRTbell molecule and combines them
using a statistical model to produce one highly accurate consensus sequence,
also called HiFi read, with base quality values. This tool powers the Circular
Consensus Sequencing workflow in SMRT Link.

### File search patterns

```yaml
ccs/v4:
  contents: ZMWs generating CCS
  max_filesize: 1024
  num_lines: 2
ccs/v5:
  contents: '"id": "ccs_processing"'
  fn: "*.json"
```
