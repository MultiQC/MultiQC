# [<img src="multiqc/templates/default/assets/img/MultiQC_logo.png" width="300" title="MultiQC">](https://github.com/ewels/MultiQC)

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

MultiQC is written in Python (v2.7 / v3.4) and contains modules for a number of common tools.
Currently, these include:

* Quality control & pre-processing
  * [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [FastQ Screen](http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
  * [Cutadapt](https://code.google.com/p/cutadapt/)
* Read aligners
  * [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
  * [STAR](https://github.com/alexdobin/STAR)
  * [Tophat](https://ccb.jhu.edu/software/tophat/)
  * [Bowtie](http://bowtie-bio.sourceforge.net)
  * [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/)
* Post-alignment tools
  * [Subread featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
  * [Picard MarkDuplicates](http://broadinstitute.github.io/picard/)
  * [Preseq](http://smithlabresearch.org/software/preseq/)
  * [Qualimap](http://qualimap.bioinfo.cipf.es/)

More to come soon. Please suggest any ideas as a new
[issue](https://github.com/ewels/MultiQC/issues) _(include an example log
file if possible)_.

## Graphical Usage

MultiQC comes with a graphical app for OS X. To use, download `MultiQC.app.zip`
from the [releases page](https://github.com/ewels/MultiQC/releases)
and unzip the archive. Double click MultiQC.app to launch, then
drag your analysis directory onto the window.

The app can be run from anywhere, though we recommend copying to your
Applications directory.

A similar graphical utility for Windows is planned for a future release.

## Command Line Usage

You can install MultiQC from [PyPI](https://pypi.python.org/pypi/multiqc/)
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
#### [v0.3.1](https://github.com/ewels/MultiQC/releases/tag/v0.3.1) - 2015-11-04
* Hotfix patch to fix broken FastQC module (wasn't finding `.zip` files properly)
* General Stats table colours now flat. Should improve browser speed.
* Empty rows now hidden if appear due to column removal in general stats
* FastQC Kmer plot removed until we have something better to show.

#### [v0.3](https://github.com/ewels/MultiQC/releases/tag/v0.3) - 2015-11-04
* Lots of lovely new documentation!
* Child templates - easily customise specific parts of the default report template
* Plugin hooks - allow other tools to execute custom code during MultiQC execution
* New Preseq module
* New design for general statistics table (snazzy new background bars)
* Further development of toolbox
  * New button to clear all filters
  * Warnings when samples are hidden, plus empty plots and table cols are hidden
  * Active toolbar tab buttons are highlighted
* Lots of refactoring by @moonso to please the Pythonic gods
  * Switched to click instead of argparse to handle command line arguments
  * Code generally conforms to best practices better now.
* Now able to supply multiple directories to search for reports
* Logging output improved (now controlled by `-q` and `-v` for quiet and verbose)
* More HTML output dealt with by the base module, less left to the modules
  * Module introduction text
  * General statistics table now much easier to add to (new helper functions)
* Images, CSS and Javascript now included in HTML, meaning that there is a single
  report file to make sharing easier
* More accessible scrolling in the report - styled scrollbars and 'to top' button.
* Modules and templates now use setuptools entry points, facilitating plugins
  by other packages. Allows niche extensions whilst keeping the core codebase clean.
* The general stats table now has a sticky header row when scrolling, thanks to
  some new javascript wizardry...
* General stats columns can have a _shared key_ which allows common colour schemes
  and data ranges. For instance, all columns describing a read count will now share
  their scale across modules.
* General stats columns can be hidden and reordered with a new modal window.
* Plotting code refactored, reports with many samples (>50 by default) don't
  automatically render to avoid freezing the browser.
* Plots with highlighted and renamed samples now honour this when exporting to
  different file types.

#### [v0.2](https://github.com/ewels/MultiQC/releases/tag/v0.2) - 2015-09-18
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

