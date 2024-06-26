---
name: MALT
url: https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/malt/
description: >
  MEGAN alignment tool
---

MALT performs alignment of metagenomic reads against a database of reference sequences (such as NR, GenBank or Silva) and produces a MEGAN RMA file as output.

The MALT MultiQC module reads the header of the MALT log files
and procudes three MultiQC sections:

- A MALT summary statistics table
- A Mappability bargraph
- A Taxonomic assignment success bargraph
