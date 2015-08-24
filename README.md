# MultiQC
MultiQC is a tool to aggregate bioinformatics results across many samples into
a single report.

It is written in Python and contains modules for a number of common tools.
Currently, these include:
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FastQ Screen](http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
* [Cutadapt](https://code.google.com/p/cutadapt/)
* [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* [STAR](https://github.com/alexdobin/STAR)
* [Subread featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
* [Picard MarkDuplicates](http://broadinstitute.github.io/picard/)

More to come soon - please suggest any ideas as a new
[issue](https://github.com/ewels/MultiQC/issues)).

## Installation

Either [download](https://github.com/ewels/MultiQC/archive/master.zip) or clone
this repository to your computer:
```
git clone https://github.com/ewels/MultiQC.git
```

Install the MultiQC python package:
```
cd MultiQC
python setup.py develop
```

## Usage

In the terminal, go to your analysis directory and run `multiqc .` - that's it!
```
cd my_analysis
multiqc .
```

The report is created in `multiqc_report/multiqc_report.html` by default.
A zip file of the report is also generated to facilitate easy transfer and sharing.

For more detailed instructions, run `multiqc -h`:

```
$ multiqc -h
usage: multiqc [-h] [-i TITLE] [-t TEMPLATE] [-o OUTPUT_DIR]
               [-m [MODULE [MODULE ...]]] [-f] [-l LEVEL] [-v]
               <analysis directory>

MultiQC is a tool to create an aggregate report summarising the results of
bioinformatics analyses across numerous samples. To run, supply with a
directory to scan for analysis results. To run here, use 'multiqc .'

positional arguments:
  <analysis directory>  Directory with analysis results to parse. Use '.' for
                        current working directory.

optional arguments:
  -h, --help            show this help message and exit
  -i TITLE, --title TITLE
                        Report title
  -t TEMPLATE, --template TEMPLATE
                        Report template to use. Choose from: default
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory. Default: multiqc_report
  -m [MODULE [MODULE ...]], --module [MODULE [MODULE ...]]
                        Use only these modules. Choose from: featureCounts,
                        bismark, star, cutadapt, fastq_screen, fastqc, picard
  -f, --force           Overwrite any existing reports
  -l LEVEL, --logging LEVEL
                        Level of logging to be printed to the console. Choose
                        from 'debug', 'info', 'warning', 'error', 'critical'
  -v, --version         Print the version of the program and exit
```

## Contributions & Support
Contributions and suggestions for new features are welcome, as are bug reports!
Please create a new [issue](https://github.com/ewels/MultiQC/issues) for any
of these, including example reports where possible.

Pull requests with new code are always gladly received, see the
[contributing notes](https://github.com/ewels/MultiQC/blob/master/CONTRIBUTING.md)
for details.

If in doubt, feel free to get in touch with the author:
[@ewels](https://github.com/ewels) (phil.ewels@scilifelab.se)
