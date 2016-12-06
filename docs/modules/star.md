---
Name: STAR
URL: https://github.com/alexdobin/STAR
Description: >
    STAR is an ultrafast universal RNA-seq aligner.
---

[STAR](https://github.com/alexdobin/STAR) is an ultrafast universal RNA-seq aligner.

This MultiQC module parses summary statistics from the `Log.final.out` log files.
Sample names are taken either from the filename prefix (`sampleNameLog.final.out`)
when set with `--outFileNamePrefix` in STAR. If there is no filename prefix,
the sample name is set as the name of the directory containing the file.

In addition to this summary log file, the module parses `ReadsPerGene.out.tab`
files generated with `--quantMode GeneCounts`, if found. 