---
name: JCVI Genome Annotation
url: https://github.com/tanghaibao/jcvi
description: A tool to compute statistics on genome annotation.
---

JCVI Genome Annotation provides quantitative measures for the assessment of genome
annotation.

The JCVI module parses the output of `python -m jcvi.annotation.stats genestats <input.gff>`.
The file name is used as the sample name.
If the output from the `python -m jcvi.annotation.stats stats <input.gff>` is present in the same directory,
it is used to draw complementary plots.

A typical result directory will contain:

```
├── Exon_Count
│   ├── sample1.txt
│   └── sample2.txt
├── Exon_Length
│   ├── sample1.txt
│   └── sample2.txt
├── Gene_Length
│   ├── sample1.txt
│   └── sample2.txt
├── Intron_Length
│   ├── sample1.txt
│   └── sample2.txt
├── sample1_genestats.txt
└── sample2_genestats.txt

```

The JCVI module has been tested with output from JCVI v1.0.9.
