---
name: vg
url: https://github.com/vgteam/vg
description: >
  VG is a suite of tools designed to allow users to manipulate graphical genomes, including aligning reads to graphical genomes.
---

The vg module parses [vg stats](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe#evaluating-with-vg-stats) alignment reports. It is specifically designed to parse the output of a vg stats report summarizing the alignment of reads to a graphical genome in a GAM file.
vg stats is capable of producing many reports summarizing many aspects of graphical genomes, including specific aspects of aligned gam files such as node coverage. This module is not meant to gather those data. Rather, this module is designed to summarize the alignment performance of gam files produced by vg giraffe created from the stdout of the vg stats command:

```sh
vg stats -a mapped.gam
```

It is not guaranteed that output created using any other parameter combination can be parsed using this module.
An example output file is below:

```txt
---
Total alignments: 727413268
Total primary: 727413268
Total secondary: 0
Total aligned: 717826332
Total perfect: 375143620
Total gapless (softclips allowed): 714388968
Total paired: 727413268
Total properly paired: 715400510
Alignment score: mean 129.012, median 132, stdev 31.5973, max 161 (244205781 reads)
Mapping quality: mean 52.8552, median 60, stdev 17.7742, max 60 (589259353 reads)
Insertions: 3901467 bp in 1466045 read events
Deletions: 6759252 bp in 2795331 read events
Substitutions: 281648245 bp in 281648245 read events
Softclips: 11480269152 bp in 252773804 read events
Total time: 291465 seconds
Speed: 2495.71 reads/second
```

For the vg module to discover the [vg stats](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe#evaluating-with-vg-stats) reports, the file must match one of the following patterns:

- "\*.txt"
- "\*.vgstats"

The graphical reports are designed to mimic a samtools stats report, including:
  1) A bar chart of the number of reads aligned/not aligned
  2) A bee swarm of:
     a) Total sequences
     b) Total properly paired alignments
     c) Total gapless alignments
     d) Total perfect alignments
     e) Mapping quality
     f) Insertions (reads)
     g) Deletions (reads)
     h) Substitutions (reads)
     i) Softclips (reads)
     j) Insertions (bases)
     k) Deletions (bp)
     l) Substitutions (bp)
     m) Softclips (bp)
     n) % of total reads aligned
     o) % of total reads properly paired

