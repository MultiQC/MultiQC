---
Name: Supernova
URL: https://www.10xgenomics.com/
Description: >
    Supernova is a de novo genome assembler for 10X Genomics linked-reads.
---

The Supernova module parses the reports from an assembly run. As a bare minimum it requires the file `report.txt`, found in the folder `sampleID/outs/`, to function. Note! If you are anything like the author (@remiolsen), you might only have files (often renamed to, e.g. `sampleID-report.txt`) lying around due to disk space limitations and for ease of sharing with your colleagues. This module will search for `*report*.txt`. If available the stats in the report file will be superseded by the higher precision numbers found in the file `sampleID/outs/assembly/stats/summary.json`. In the same folder, this module will search for the following plots and render them:

* `histogram_molecules.json` -- Inferred molecule lengths
* `histogram_kmer_count.json` -- Kmer multiplicity

This module has been tested for Supernova versions `1.1.4` and `1.2.0`
