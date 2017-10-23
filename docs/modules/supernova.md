---
Name: Supernova
URL: https://www.10xgenomics.com/
Description: >
    Supernova is a de novo genome assembler for 10X Genomics linked-reads.
---

#### Important notes

Due to the size of the `histogram_kmer_count.json` files, MultiQC is likely to skip these files. To be able to display these you will need to change the MultiQC configuration to allow for larger logfiles, see the MultiQC [documentation](http://multiqc.info/docs/#troubleshooting-1). For instance, if you run MultiQC as part of an analysis pipeline, you can create a `multiqc_config.yaml` file in the working directory, containing the following line:

```
log_filesize_limit: 100000000
```

#### General Notes

The Supernova module parses the reports from an assembly run. As a bare minimum it requires the file `report.txt`, found in the folder `sampleID/outs/`, to function. Note! If you are anything like the author (@remiolsen), you might only have files (often renamed to, e.g. `sampleID-report.txt`) lying around due to disk space limitations and for ease of sharing with your colleagues. This module will search for `*report*.txt`. If available the stats in the report file will be superseded by the higher precision numbers found in the file `sampleID/outs/assembly/stats/summary.json`. In the same folder, this module will search for the following plots and render them:

* `histogram_molecules.json` -- Inferred molecule lengths
* `histogram_kmer_count.json` -- Kmer multiplicity

This module has been tested using Supernova versions `1.1.4` and `1.2.0`

