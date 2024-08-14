---
name: BioBloom Tools
urls: ['https://github.com/bcgsc/biobloom/']
summary: >
  Assigns reads to different references using bloom filters. This is faster than alignment and can be used for contamination detection
extra_description: >
  BioBloom tools (BBT) create filters for a given reference and then to categorize sequences.
  This methodology is faster than alignment but does not provide mapping locations. BBT was initially intended to
  be used for pre-processing and QC applications like contamination detection, but is flexible to accommodate other
  purposes. This tool is intended to be a pipeline component to replace costly alignment steps.
---

### File search patterns

```yaml
biobloomtools:
  contents: "filter_id\thits\tmisses\tshared\trate_hit\trate_miss\trate_shared"
  num_lines: 2
```