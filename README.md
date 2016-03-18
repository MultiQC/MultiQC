<img src="multiqc/templates/default/assets/img/MultiQC_logo.png" width="300" title="MultiQC">

### Aggregate bioinformatics results across many samples into a single report.

##### Find [documentation](http://multiqc.info/docs/0.2/README.md) and [example reports](http://multiqc.info/examples/rna-seq/multiqc_report.html) at [multiqc.info](http://multiqc.info)

<table>
  <tr>
    <th>Stable:</th>
    <td>
      <a title="PyPI Version" href="https://pypi.python.org/pypi/multiqc/"><img src="https://img.shields.io/pypi/v/multiqc.svg"></a>
      <a title="Stable docs" href="http://multiqc.info/docs/"><img src="https://img.shields.io/badge/docs-stable-green.svg"></a>
      <a title="Licence" href="https://github.com/ewels/MultiQC/blob/master/LICENSE"><img src="https://img.shields.io/pypi/l/multiqc.svg"></a>
      <a title="PyPI Downloads" href="http://multiqc.info/stats.php"><img src="https://img.shields.io/pypi/dm/multiqc.svg"></a>
    </td>
  </tr>
  <tr>
    <th>Devel:</th>
    <td>
      <a title="Build Status" href="https://travis-ci.org/ewels/MultiQC"><img src="https://travis-ci.org/ewels/MultiQC.svg?branch=master"></a>
      <a title="Devel docs" href="https://github.com/ewels/MultiQC/tree/master/docs"><img src="https://img.shields.io/badge/docs-devel-yellow.svg"></a>
      <img src="https://img.shields.io/badge/Python-2.7-green.svg">
      <img src="https://img.shields.io/badge/Python-3.4-green.svg">
    </td>
  </tr>
</table>

-----

MultiQC is a tool to create a single report with interactive plots
for multiple bioinformatics analyses across many samples.

MultiQC is written in Python (tested with v2.7 / v3.4 / v3.5) and is
available on the [Python Package Index](https://pypi.python.org/pypi/multiqc/).

Currently, supported tools include:

* Quality control & pre-processing
  * [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [FastQ Screen](http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
  * [Cutadapt](https://code.google.com/p/cutadapt/)
  * [Skewer](https://github.com/relipmoc/skewer)
* Read aligners
  * [STAR](https://github.com/alexdobin/STAR)
  * [Tophat](https://ccb.jhu.edu/software/tophat/)
  * [Bowtie](http://bowtie-bio.sourceforge.net)
  * [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/)
  * [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
  * [HiCUP](http://www.bioinformatics.babraham.ac.uk/projects/hicup/)
* Post-alignment tools
  * [Subread featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
  * [Samblaster](https://github.com/GregoryFaust/samblaster)
  * [Picard](http://broadinstitute.github.io/picard/) (MarkDuplicates, InsertSizeMetrics, GcBiasMetrics, HsMetrics, OxoGMetrics)
  * [Preseq](http://smithlabresearch.org/software/preseq/)
  * [Qualimap](http://qualimap.bioinfo.cipf.es/) (BamQC, RNASeq)

More modules are being written all of the time. Please suggest any ideas as a new
[issue](https://github.com/ewels/MultiQC/issues) _(include an example log
file if possible)_.

## Installation

You can install MultiQC from [PyPI](https://pypi.python.org/pypi/multiqc/)
using `pip` as follows:

```
pip install multiqc
```

If you would like the development version instead, the command is:

```
pip install git+https://github.com/ewels/MultiQC.git
```

Alternatively, you can install using [Conda](http://anaconda.org/)
from the [bioconda channel](https://bioconda.github.io/):
```
conda install -c bioconda multiqc
```

## Usage
Once installed, you can use MultiQC by navigating to your analysis directory
(or a parent directory) and running the tool:

```
multiqc .
```

That's it! MultiQC will scan the specified directory ('.' is the current dir)
and produce a report detailing whatever it finds.

The report is created in `multiqc_report.html` by default. Tab-delimited data
files are also created in `multiqc_data/`, containing extra information.
These can be easily inspected using Excel (use `--data-format` to get `yaml`
or `json` instead).

For more detailed instructions, run `multiqc -h` or see the
[documentation](http://multiqc.info/docs/#running-multiqc).

## Contributions & Support

Contributions and suggestions for new features are welcome, as are bug reports!
Please create a new [issue](https://github.com/ewels/MultiQC/issues) for any
of these, including example reports where possible.

Pull requests with new code are always gladly received, see the
[contributing notes](https://github.com/ewels/MultiQC/blob/master/CONTRIBUTING.md)
for details. These notes include extensive help with how to use the built in code.

If in doubt, feel free to get in touch with the author:
[@ewels](https://github.com/ewels) (phil.ewels@scilifelab.se)

### Contributors
Many thanks to those who have helped out with with project.
* Project lead and main author: @ewels
* Early code refactoring help: @moonso
* Early version of Qualimap module: @guillermo-carrasco
* Skewer and Samblaster modules: @dakl
* Samtools stats module: @lpantano
* Tweaks / core help: @robinandeer and @avilella

