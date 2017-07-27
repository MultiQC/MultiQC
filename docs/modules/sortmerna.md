---
Name: SortMeRNA
URL: http://bioinfo.lifl.fr/RNA/sortmerna/
Description: >
    SortMeRNA is a program tool for filtering, mapping and OTU-picking NGS reads
    in metatranscriptomic and metagenomic data.
---

SortMeRNA is a program tool for filtering, mapping and OTU-picking NGS reads in metatranscriptomic and metagenomic data. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of nucleotide sequences. The main application of SortMeRNA is filtering ribosomal RNA from metatranscriptomic data.

The MultiQC module parses the log files, which are created when `SortMeRNA` is run with the `--log` option.

The default header in the 'General Statistics' table is '% rRNA'. Users can override this using the configuration option:

```
sortmerna:
    tab_header: 'My database hits'
```