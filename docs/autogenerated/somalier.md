---
name: Somalier
urls: ["https://github.com/brentp/somalier"]
summary: >
  Genotype to pedigree correspondence checks from sketches derived from BAM/CRAM or VCF
---

Somalier can be used to find sample swaps or duplicates in cancer
projects, where there is often no jointly-called VCF across samples.
It is also extremely efficient and so can be used to find all-vs-all
relatedness estimates for thousands of samples.
It also outputs information on sex, depth, heterozgyosity, and ancestry
to be used for general QC.

### File search patterns

```yaml
somalier/pairs:
  contents: hom_concordance
  fn: "*.pairs.tsv"
  num_lines: 5
somalier/samples:
  contents: "#family_id"
  fn: "*.samples.tsv"
  num_lines: 5
somalier/somalier-ancestry:
  fn: "*.somalier-ancestry.tsv"
```
