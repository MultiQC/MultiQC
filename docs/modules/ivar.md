---
Name: iVar
URL: https://github.com/andersen-lab/ivar
Description: >
    Functions for viral amplicon-based sequencing.
---

iVar is a computational package that contains functions broadly useful for viral amplicon-based sequencing.

This module parses the output from the `ivar trim` command and creates a table view. Both output from V1 and V2 of the tool are supported and parsed accordingly. The executable used can easily be installed from the BioConda channel using `conda install -c bioconda ivar=1.2` for example. The barplot shows a simple classification of reads that are passing criteria, are too short for further usage or outside of primer regions.
