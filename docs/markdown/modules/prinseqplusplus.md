---
name: PRINSEQ++
url: https://github.com/Adrian-Cantu/PRINSEQ-plus-plus
description: PRINSEQ++ is a C++ implementation of the prinseq-lite.pl program.
---

PRINSEQ++ can be used to filter, reformat or trim genomic and metagenomic sequence data. It is 5X faster than `prinseq-lite.pl` and uses less RAM thanks to the use of multi-threading and the `cboost` libraries. It can read and write compressed (gzip) files, drastically reducing the use of hard drive.

This module requires that PRINSEQ++ has been run with the flag `-VERBOSE 1`.
It uses the log file name as the sample name.
