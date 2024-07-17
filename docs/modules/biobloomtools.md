---
name: BioBloom Tools
url: https://github.com/bcgsc/biobloom/
description: Assigns reads to different references using bloom filters. This is faster than alignment and can be used for contamination detection
---

BioBloom tools (BBT) create filters for a given reference and then to categorize sequences.
This methodology is faster than alignment but does not provide mapping locations. BBT was initially intended to
be used for pre-processing and QC applications like contamination detection, but is flexible to accommodate other
purposes. This tool is intended to be a pipeline component to replace costly alignment steps.
