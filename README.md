# [<img src="multiqc/templates/default/assets/img/MultiQC_logo.png" width="300" title="MultiQC">](https://github.com/ewels/MultiQC)

**MultiQC is a tool to aggregate bioinformatics results across many samples into
a single report.**

[![Build Status](https://travis-ci.org/ewels/MultiQC.svg?branch=master)](https://travis-ci.org/ewels/MultiQC)

See an example report [here](http://multiqc.info/examples/rna-seq/multiqc_report.html) (more at [multiqc.info](http://multiqc.info)).

It is written in Python and contains modules for a number of common tools.
Currently, these include:

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FastQ Screen](http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
* [Cutadapt](https://code.google.com/p/cutadapt/)
* [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* [STAR](https://github.com/alexdobin/STAR)
* [Tophat](https://ccb.jhu.edu/software/tophat/)
* [Bowtie](http://bowtie-bio.sourceforge.net)
* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/)
* [Subread featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
* [Picard MarkDuplicates](http://broadinstitute.github.io/picard/)

More to come soon. Please suggest any ideas as a new
[issue](https://github.com/ewels/MultiQC/issues).

## Graphical Usage

MultiQC comes with a graphical app for OS X. To use, download `MultiQC.app.zip`
from the [releases page](https://github.com/ewels/MultiQC/releases)
and unzip the archive. Double click MultiQC.app to launch, then
drag your analysis directory onto the window.

The app can be run from anywhere, though we recommend copying to your
Applications directory.

A similar graphical utility for Windows is planned for a future release.

## Command Line Usage

You can install MultiQC from [PyPI](https://pypi.python.org/pypi/multiqc/0.1)
using `pip` as follows:

```
pip install multiqc
```

If you would like the development version instead, the command is:

```
pip install git+https://github.com/ewels/MultiQC.git
```

Then it's just a case of going to your analysis directory and running the script:

```
multiqc .
```

That's it! MultiQC will scan the specified directory ('.' is the current dir)
and produce a report detailing whatever it finds.

The report is created in `multiqc_report/multiqc_report.html` by default.
A zip file of the report is also generated to facilitate easy transfer and sharing.

Tab-delimited data files are also created in `multiqc_report/report_data/`,
containing extra information. These can be easily inspected using Excel.

For more detailed instructions, run `multiqc -h`

## Contributions & Support

Contributions and suggestions for new features are welcome, as are bug reports!
Please create a new [issue](https://github.com/ewels/MultiQC/issues) for any
of these, including example reports where possible.

Pull requests with new code are always gladly received, see the
[contributing notes](https://github.com/ewels/MultiQC/blob/master/CONTRIBUTING.md)
for details. These notes include extensive help with how to use the built in code.

If in doubt, feel free to get in touch with the author:
[@ewels](https://github.com/ewels) (phil.ewels@scilifelab.se)

## Version History
#### v0.2 - 2015-09-18
* Code restructuring for nearly all modules. Common base module
  functions now handle many more functions (plots, config, file import)
  * See the [contributing notes](https://github.com/ewels/MultiQC/blob/master/CONTRIBUTING.md)
    for instructions on how to use these new helpers to make your own module
* New report toolbox - sample highlighting, renaming, hiding
  * Config is autosaved by default, can also export to a file for sharing
  * Interactive tour to help users find their way around
* New Tophat, Bowtie 2 and QualiMap modules
  * Thanks to @guillermo-carrasco for the QualiMap module
* Bowtie module now works
* New command line parameter `-d` prefixes sample names with the directory that
  they were found in. Allows duplicate filenames without being overwritten.
* Introduction walkthrough helps show what can be done in the report
* Now compatible with both Python 2 and Python 3
* Software version number now printed on command line properly, and in reports.
* Bugfix: FastQC doesn't break when only one report found
* Bugfix: FastQC seq content heatmap highlighting
* Many, many small bugfixes

#### [v0.1](https://github.com/ewels/MultiQC/releases/tag/v0.1) - 2015-09-01
* The first public release of MultiQC, after a month of development. Basic
structure in place and modules for FastQC, FastQ Screen, Cutadapt, Bismark, 
STAR, Bowtie, Subread featureCounts and Picard MarkDuplicates. Approaching
stability, though still under fairly heavy development.

