---
Name: Peddy
URL: https://github.com/brentp/peddy
Description: >
    Peddy calculates genotype :: pedigree correspondence checks,
    ancestry checks and sex checks using VCF files.
---

Peddy compares familial-relationships and sexes as reported in a
PED file with those inferred from a VCF.

It samples the VCF at about 25000 sites (plus chrX) to accurately
estimate relatedness, IBS0, heterozygosity, sex and ancestry.
It uses 2504 thousand genome samples as backgrounds to calibrate
the relatedness calculation and to make ancestry predictions.

It does this very quickly by sampling, by using C for computationally
intensive parts, and by parallelization.