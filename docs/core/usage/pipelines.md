---
title: Using MultiQC in pipelines
description: Integration within your workflow manager of choice
---

# Using MultiQC in pipelines

MultiQC has been designed to be placed at the end of bioinformatics workflows
and works well when it's a final step in a pipeline. Here you can find a few
tips about integration with different workflow tools.

I'll use FastQC as an example input in all of these examples as it's a common
use case. However, the concepts should apply for any of the tools that MultiQC supports.

:::tip
Remember that you can use [Custom Content](../custom_content/index.md)
feature to easily collect pipeline-specific metadata (pipeline run-time data,
links to documentation, etc.) in to a format that can be inserted in to your report.
:::

:::tip
Want to include software version information in your report?
Even group versions by the step in the pipeline? You can! 🎉
See the [Listing software versions](../reports/customisation.md#listing-software-versions) documentation.
:::

If you know exactly which modules will be used by MultiQC, you can use the
`-m`/`--modules` flag to specify just these. This will speed up MultiQC a little.
This will probably only make a noticeable impact if your pipeline has thousands
of input files for MultiQC.

In the following, we show examples for [Nextflow](#nextflow) and [Snakemake](#snakemake) integration.

## Nextflow

An example mini-pipeline for [nextflow](https://www.nextflow.io/) which runs FastQC
and MultiQC is below.
See the [nf-core](https://nf-co.re/) pipelines for lots of examples of full nextflow
pipelines that use MultiQC.

```groovy
#!/usr/bin/env nextflow

params.reads = "data/*{1,2}.fastq.gz"
Channel.fromFilePairs( params.reads ).into { read_files_fastqc }

process fastqc {
    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

process multiqc {
    input:
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"

    script:
    """
    multiqc .
    """
}
```

You will need to specify channels for all files that MultiQC needs to run, so that nextflow
stages them correctly in the MultiQC work directory.

Note that `.collect()` is needed to make MultiQC run once for all upstream outputs.

The `.ifEmpty([])` add on isn't really needed here, but is helpful in larger pipelines where
some processes may be optional. Without this, if any channels are not processed then MultiQC
won't run.

### Clashing input filenames

If a nextflow process tries to stage more than one input file with an identical filename,
it will throw an error. Putting inputs into their own subfolder (`file ('fastqc/*')`) in
the above examples is not really needed, but reduces the chance of input filename clashes.

If you're using a tool that gives the same filename to each file that MultiQC uses, you'll
need to tell nextflow to rename the inputs to prevent clashes.

For example, [StringTie](https://ccb.jhu.edu/software/stringtie/) prints statistics to
STDOUT that MultiQC uses to generate reports. We can easily collect this in nextflow by
using the `.command.log` file that nextflow saves as an output file from the process
(`file ".command.log" into stringtie_log`). However, now every sample has the same filename
for MultiQC.

We get around this by using [dynamic input file names](https://www.nextflow.io/docs/latest/process.html#dynamic-input-file-names)
with nextflow:

```groovy
file ('stringtie/stringtie_log*') from stringtie_log.collect().ifEmpty([])
```

This `file` pattern renames each stringtie log file to `stringtie_log1`,
`stringtie_log2` and so on, ensuring that we avoid any filename clashes.

Note that MultiQC finds output from some tools based on their filename, so use with caution
(you may need to define some custom [module search patterns](../getting_started/config.md#module-search-patterns)).

### Custom run name

When you launch nextflow, you can use the `-name` command line flag (single hyphen) to give
a name to that specific pipeline run. You can configure your pipeline to pass this on to
MultiQC. It can then be used as the report title and filename.

Here's an example snippet for this:

```groovy
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

process multiqc {
    input:
    file ('fastqc/*') from ch_fastqc_results_for_multiqc.collect().ifEmpty([])

    output:
    file "*_report.html" into ch_multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc $rtitle $rfilename .
    """
}
```

Note that we use `params.name` as a placeholder, this gives the benefit that both `-name`
and the common typo `--name` work for this case.

There isn't an easy way to know if a custom value for `-name` has been given to nextflow,
but all default names are two lowercase words with a single underscore so if the name
matches this pattern then we ignore it.

### MultiQC config file

It can be nice to use a config file for MultiQC to add in some static content to reports
about your pipeline (eg. `report_comment`).
Remember that even this config file should also be in a nextflow channel,
so that nextflow correctly stages it (especially important when running on the cloud).

This snippet works with a `params` variable again, so that pipeline users can replace
the config file with one of their own if they wish.

```groovy
params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"
Channel.fromPath(params.multiqc_config, checkIfExists: true).set { ch_config_for_multiqc }

process multiqc {
    input:
    file multiqc_config from ch_config_for_multiqc
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"

    script:
    """
    multiqc --config $multiqc_config .
    """
}
```

## Snakemake

The best-practice for using MultiQC in [Snakemake](https://snakemake.readthedocs.io/) pipelines is to use a predefined [Snakemake wrapper](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/multiqc.html). For example, see the following mini-pipeline:

```yaml
configfile: "config.yaml"

rule all:
  input: "multiqc_report.html"

rule fastqc:
  input: "data/{sample}.fastq.gz"
  output: html="fastqc/{sample}.html",
    zip="fastqc/{sample}.zip"
  wrapper: "0.31.1/bio/fastqc"

rule multiqc:
  input: expand("fastqc/{sample}.html", sample=config["samples"])
  output: "multiqc_report.html"
  wrapper: "0.31.1/bio/multiqc"
```

Snakemake wrappers not only deliver predefined and unit tested code for generating the requested output with the respective tool, but also define the required software stack in terms of conda packages.

You can run the above workflow as follows:

```bash
snakemake --use-conda
```

This first installs all the required tools into isolated conda environments, and then executes all necessary steps to create the target that is given in the top rule. In other words, it will execute FastQC twice (creating `fastqc/1.html` and `fastqc/2.html`) and MultiQC once, creating `multiqc_report.html`.

Of course, using conda is optional, but it greatly increases reproducibility.
Snakemake is not limited to wrappers (although its [wrapper repository](https://snakemake-wrappers.readthedocs.io) provides many in the field of bioinformatics), but also supports direct execution of [shell commands](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rules) and integration of [custom scripts](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts) (e.g., for plotting).

If you prefer to use MultiQC without a snakemake wrapper, you can see a minimal example on GitHub: [jakevc/snakemake_multiqc](https://github.com/jakevc/snakemake_multiqc). This has an example script and some test data for you to play with.
