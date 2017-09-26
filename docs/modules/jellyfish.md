---
Name: Jellyfish
URL: http://www.cbcb.umd.edu/software/jellyfish/
Description: >
    JELLYFISH is a tool for fast, memory-efficient counting of k-mers in DNA.
---

JELLYFISH is a tool for fast, memory-efficient counting of k-mers in DNA. A k-mer is a substring of length k, and counting the occurrences of all such substrings is a central step in many analyses of DNA sequence. JELLYFISH can count k-mers using an order of magnitude less memory and an order of magnitude faster than other k-mer counting packages by using an efficient encoding of a hash table and by exploiting the "compare-and-swap" CPU instruction to increase parallelism.

The MultiQC module for Jellyfish parses *only* `*_jf.hist` files. The general usage of jellyfish to be parsed by MultiQC module needs to be:

 - `gunzip -c file.fastq.gz | jellyfish count -o file.jf  -m ...`
 - `jellyfish histo -o file_jf.hist -f file.jf`

In case a user wants to customise the matching pattern for jellyfish, then multiqc can be run with the option `--cl_config "sp: { jellyfish: { fn: 'PATTERN' } }"` where `PATTERN` is the pattern to be matched. For example:

```
multiqc . --cl_config "sp: { jellyfish: { fn: '*.hist' } }"
```
