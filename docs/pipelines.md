# Using MultiQC in pipelines

MultiQC has been designed to be placed at the end of bioinformatics workflows
and works well when it's a final step in a pipeline. Here you can find a few
tips about integration with different workflow tools.

I'll use FastQC as an example input in all of these examples as it's a common
use case. However, the concepts should apply for any of the tools that MultiQC supports.

Remember that you can use [Custom Content](https://multiqc.info/docs/#custom-content)
feature to easily collect pipeline-specific metadata (software version numbers,
pipeline run-time data, links to documentation) in to a format that can be inserted
in to your report.

If you know exactly which modules will be used by MultiQC, you can use the
`-m`/`--modules` flag to specify just these. This will speed up MultiQC a little.
This will probably only make a noticeable impact if your pipeline has thousands
of input files for MultiQC.

## Nextflow

An example mini-pipeline for [nextflow](https://www.nextflow.io/) which runs FastQC
and MultiQC is below.
See the [nf-core](https://nf-co.re/) pipelines for lots of examples of full nextflow
pipelines that use MultiQC.

```groovy
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
and the common type `--name` work for this case.

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
