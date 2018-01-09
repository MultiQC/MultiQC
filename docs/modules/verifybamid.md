---
Name: verify BAM ID
URL: https://genome.sph.umich.edu/wiki/VerifyBamID
Description: >
verifyBamID is a software that verifies whether the reads in particular file match previously known genotypes for an individual (or group of individuals), and checks whether the reads are contaminated as a mixture of two samples. 
verifyBamID can detect sample contamination and swaps when external genotypes are available. When external genotypes are not available, verifyBamID still robustly detects sample swaps.
---
A key step in any genetic analysis is to verify whether data being generated matches expectations. verifyBamID checks whether reads in a BAM file match previous genotypes for a specific sample. In addition, it detects possible sample mixture from population allele frequency only, which can be particularly useful when the genotype data is not available.

Using a mathematical model that relates observed sequence reads to an hypothetical true genotype, verifyBamID tries to decide whether sequence reads match a particular individual or are more likely to be contaminated (including a small proportion of foreign DNA), derived from a closely related individual, or derived from a completely different individual.