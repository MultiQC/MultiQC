---
name: CheckQC
urls: ['https://github.com/Molmed/checkQC']
summary: >
  Checks a set of quality criteria against an Illumina runfolder
---

The module parses a CheckQC JSON file, so make sure to use CheckQC with the `--json` flag and collect the stdout in a file.

### File search patterns

```yaml
checkqc:
  contents: instrument_and_reagent_type
  fn: '*.json'
```