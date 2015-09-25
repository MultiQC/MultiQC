---
title: Installation
layout: toc
---

# Installing MultiQC

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


## Graphical Usage

MultiQC comes with a graphical app for OS X. To use, download `MultiQC.app.zip`
from the [releases page](https://github.com/ewels/MultiQC/releases)
and unzip the archive. Double click MultiQC.app to launch, then
drag your analysis directory onto the window.

The app can be run from anywhere, though we recommend copying to your
Applications directory.

A similar graphical utility for Windows is planned for a future release.