---
Name: GATK
URL: https://www.broadinstitute.org/gatk/
Description: Variant Discovery in High-Throughput Sequencing Data
---

Developed by the [Data Science and Data Engineering](http://www.broadinstitute.org/dsde)
group at the [Broad Institute](http://www.broadinstitute.org/), the GATK toolkit offers
a wide variety of tools with a primary focus on variant discovery and genotyping.
Its powerful processing engine and high-performance computing features make it capable
of taking on projects of any size.

Supported tools: 
- `BaseRecalibrator`
- `VariantEval`

#### BaseRecalibrator
[BaseRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)
is a tool for detecting systematic errors in read base quality scores of aligned high-throughput
sequencing reads. It outputs a base quality score recalibration table that can be used in 
conjunction with the 
[PrintReads](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php)
tool to recalibrate base quality scores.

### VariantEval
[VariantEval](https://software.broadinstitute.org/gatk/gatkdocs/current/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php)
is a general-purpose tool for variant evaluation. It gives information about percentage of
variants in dbSNP, genotype concordance, Ti/Tv ratios and a lot more.
