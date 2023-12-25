---
name: biobambam2
url: https://github.com/gt1/biobambam2
description: >
  biobambam2 contains tools for processing BAM files for early stage alignment file processing
---

Currently, the biobambam2 module only processes output from the `bamsormadup` command.
Not only that, but it cheats by using the module code from Picard/MarkDuplicates.
The output is so similar that the code simply sets up a module with unique name and
filename search pattern and then uses the parsing code from the Picard module.

Apart from behind the scenes coding, this module should work in exactly the same way
as all other MultiQC modules.
