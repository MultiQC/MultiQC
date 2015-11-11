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
  * [calculate_read_overlap](http://github.com/avilella/calculate_clip_overlap/)
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

