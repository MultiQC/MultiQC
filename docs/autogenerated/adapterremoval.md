---
name: Adapter Removal
urls: ["https://github.com/mikkelschubert/adapterremoval"]
summary: >
  Removes adapter sequences, trims low quality bases from 3' ends, or merges overlapping pairs into consensus
---

AdapterRemoval searches for and removes remnant adapter sequences from High-Throughput Sequencing (HTS) data and (optionally) trims low quality bases from the 3' end of reads following adapter removal. It can analyze both single end and paired end data, and can be used to merge overlapping paired-ended reads into (longer) consensus sequences. Additionally, the AdapterRemoval may be used to recover a consensus adapter sequence for paired-ended data, for which this information is not available.

The module parses `*.settings` logs from Adapter Removal.

Supported setting file results:

- `single end`
- `paired end noncollapsed`
- `paired end collapsed`

### File search patterns

```yaml
adapterremoval:
  contents: AdapterRemoval
  fn: "*.settings"
  num_lines: 1
```