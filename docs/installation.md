# MultiQC Installation and Usage

## Installing MultiQC
You can install MultiQC from [PyPI](https://pypi.python.org/pypi/multiqc/0.1)
using `pip` (the Python package manager) as follows:
```
pip install multiqc
```

If you would like the development version instead, the command is:
```
pip install git+https://github.com/ewels/MultiQC.git
```

If you don't have `pip` installed, you can clone and install the code yourself:
```
git clone https://github.com/ewels/MultiQC.git
cd MultiQC
python setup.py install
```

## Basic Usage
Once installed, just go to your analysis directory and run `multiqc`, followed
by a list of directories to search. At it's simplest, this can just be `.`
(the current working directory):
```
multiqc .
```

That's it! MultiQC will scan the specified directories and produce a report
based on details found in any log files that it recognises.

For a description of all command line parameters, run `multiqc --help`

### Renaming reports
The report is called `multiqc_report.html` by default. Tab-delimited data files
are created in `multiqc_data/`, containing additional information.
You can prefix a custom name to these using the `-p`/`--prefix` parameter, or instruct
MultiQC to create them in a subdirectory using the `-o`/`-outdir` parameter.

### Overwriting existing reports
It's quite common to repeatedly create new reports as new analysis results
are generated. Instead of manually deleting old reports, you can just specify
the `-f` parameter and MultiQC will overwrite any conflicting report filenames.

### Sample names prefixed with directories
Sometimes, the same samples may be processed in different ways. If MultiQC
finds log files with the same sample name, the previous data will be overwritten
(this can be inspected by running MultiQC with `-v`/`--verbose`).

To avoid this, run MultiQC with the `-d`/`--dirs` parameter. This will prefix every
sample name with the directory path that the log file was found within. As
such, sample names will no longer be unique, and data will not be overwritten.


## Graphical Mac OSX App
MultiQC comes with a graphical app for OS X. See the
[MultiQC_OSXApp](https://github.com/ewels/MultiQC_OSXApp) repository
for more details and downloads.