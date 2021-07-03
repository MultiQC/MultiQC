---
Name: bclconvert
URL: https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html
Description: >
  bclconvert can be used to both demultiplex data and convert BCL files to
  FASTQ file formats for downstream analysis.
---

This BclConvert module is based on the bcl2fastq multiqc module. It can parse multiple
bclconvert run outputs as long as they are from the same sequencing run. When doing this,
the undetermined reads will be 'corrected' and re-calculated (as an unknown read from
some one bclcovnert run might not be truly unknown, but simply from another run).

#### Calculate estimated depth

You can specify a genome size in config

It's often useful to talk about sequencing yield in terms of estimated depth of coverage.
In order to make MultiQC show the estimated depth for each sample, specify the reference genome/target size in your [MultiQC configuration](http://multiqc.info/docs/#configuring-multiqc):

```yaml
bclconvert:
  genome_size: 3049315783
```

The coverage depth will be estimated as the yield Q30 dvivided by the genome size.

MultiQC comes with effective genome size presets for Human and Mouse, so you can
provide the genome build name instead, like this: `genome_size: hg38_genome`. The
following values are supported: `hg19_genome`, `hg38_genome`, `mm10_genome`.
