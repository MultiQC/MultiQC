---
title: Quick Start
description: A tutorial with a typical setup for people in a rush.
---

# MultiQC: Quick start

This tutorial covers installation and first run for a typical user.
It's not meant to be comprehensive - see the rest of the documentation for that -
it's just to get the majority up and running quickly so you can get a taste for how to use MultiQC.

## Install Conda

In order to install MultiQC, we first need Python.
Arguably, the easiest way to do this is with Conda
(see [full docs](../installation/#python-with-conda)).

1. [Download miniconda](https://conda.io/miniconda.html) for your operating system.
2. Run the bash script and follow the prompts.
3. Restart your terminal shell.
4. [Configure your conda channels](https://bioconda.github.io/#usage) to work with BioConda:
   ```bash
   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   ```
5. Create a new conda environment:
   ```bash
   conda create --name myenv python=3.11
   conda activate myenv
   ```

## Install MultiQC

Now that we have Python, we can install MultiQC.
As we're already using Conda, we may as well install MultiQC with Conda too (see [full docs](../installation/#conda)).

```bash
conda install multiqc
```

Check that it worked by printing the MultiQC version (or `--help` text):

```bash
multiqc --version
```

```txt
multiqc, version 1.20
```

## Get some example data

To try MultiQC out quickly, you can fetch some example input data from the [Example reports](https://multiqc.info/example-reports/) page.

Each example report has a link to _Download input data_.
You should be able to recreate the example report using this.

For example, for the [RNA-seq report](https://multiqc.info/example-reports/rna-seq/):

```bash
curl -O -J -L http://multiqc.info/examples/rna-seq/data.zip
unzip data.zip
```

You should now have a directory called `data` which is full of analysis result files.

For the RNA-seq example, we have outputs from a bioinformatics analysis of some publicly available data.
We have logs and reports from [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [TrimGalore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) ([Cutadapt](https://cutadapt.readthedocs.io/)), [STAR](https://github.com/alexdobin/STAR) and [featureCounts](https://subread.sourceforge.net/).

## Run MultiQC

There isn't much to running MultiQC really - just point it at the directory that contains your files and it will search recursively for anything it recognises.

Assuming that you are still in the directory where you just extracted the data,
the current working directory (`.`) contains your files:

```bash
multiqc .
```

```txt
  /// MultiQC 🔍 | v1.20

|           multiqc | Search path : /demo/data
|         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 121/121
|    feature_counts | Found 8 reports
|              star | Found 8 reports
|          cutadapt | Found 16 reports
|            fastqc | Found 32 reports
|           multiqc | Report      : multiqc_report.html
|           multiqc | Data        : multiqc_data
|           multiqc | MultiQC complete
```

## Open the report

You can see in the log output that MultiQC created a file called `multiqc_report.html`.
Open it and take a look (you can usually <kbd>ctrl</kbd>/<kbd>cmd</kbd> + click the filename in most terminals).
It should look basically the same as the example report [on the MultiQC website](https://multiqc.info/example-reports/rna-seq/).

Try using the toolbox features in the right hand sidebar, for example hiding and highlighting specific samples.

Also have a look at the directory `multiqc_data` that was created.
This contains the parsed data in a nice friendly format, ready for any further downstream analysis.
