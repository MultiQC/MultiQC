---
name: BBDuk
url: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
description: Tool for common data-quality-related trimming, filtering, and masking operations
---

The BBDuk module produces summary statistics from the stdout logging information
from the BBDuk tool of the [BBTools](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/) suite of tools.

"Duk" stands for Decontamination Using Kmers. BBDuk was developed to combine
most common data-quality-related trimming, filtering, and masking operations
into a single high-performance tool.

The module can summarise data from the following BBDuk funtionality
(descriptions from command line help output):

- `entropy` - entropy filtering
- `ktrim` - kmer trimming
- `qtrim` - quality trimming
- `maq` - read quality filtering
- `ref` contaminant filtering

Additional information on the BBMap tools is available on
[SeqAnswers](http://seqanswers.com/forums/showthread.php?t=41057).
