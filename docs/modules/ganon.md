---
name: ganon
url: https://pirovc.github.io/ganon/
description: >
  ganon is developed for, but not limited, to the metagenomics classification problem: quickly assign sequence fragments to their closest reference among thousands of references via  Interleaved Bloom Filters of k-mer/minimizers.
---

ganon is designed to index large sets of genomic reference sequences and to classify reads against them efficiently. The tool uses Interleaved Bloom Filters as indices based on k-mers/minimizers. It was mainly developed, but not limited, to the metagenomics classification problem: quickly assign sequence fragments to their closest reference among thousands of references.

The module takes summary statistics from the a file containing the re-directed stdout from `ganon classify`.
