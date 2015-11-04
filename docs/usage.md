# Basic Usage
Once installed, just go to your analysis directory and run `multiqc`, followed
by a list of directories to search. At it's simplest, this can just be `.`
(the current working directory):
```
multiqc .
```

That's it! MultiQC will scan the specified directories and produce a report
based on details found in any log files that it recognises.

See [Using MultiQC Reports](reports.md) for more information about how
to use the generated report.

For a description of all command line parameters, run `multiqc --help`.

## Renaming reports
The report is called `multiqc_report.html` by default. Tab-delimited data files
are created in `multiqc_data/`, containing additional information.
You can prefix a custom name to these using the `-p`/`--prefix` parameter, or instruct
MultiQC to create them in a subdirectory using the `-o`/`-outdir` parameter.

Note that different MultiQC templates may have different defaults.

## Overwriting existing reports
It's quite common to repeatedly create new reports as new analysis results
are generated. Instead of manually deleting old reports, you can just specify
the `-f` parameter and MultiQC will overwrite any conflicting report filenames.

## Sample names prefixed with directories
Sometimes, the same samples may be processed in different ways. If MultiQC
finds log files with the same sample name, the previous data will be overwritten
(this can be inspected by running MultiQC with `-v`/`--verbose`).

To avoid this, run MultiQC with the `-d`/`--dirs` parameter. This will prefix every
sample name with the directory path that the log file was found within. As
such, sample names will no longer be unique, and data will not be overwritten.

## Using different templates
MultiQC is built around a templating system. You can produce reports with
different styling by using the `-t`/`--template` option. The available templates
are listed with `multiqc --help`.

If you're interested in creating your own custom template, see the
[writing new templates](templates.md) section.

## Choosing which modules to run
Sometimes, it's desirable to choose which MultiQC modules run. This could be
because you're only interested in one type of output and want to keep the
reports small. Or perhaps the output from one module is misleading in your
situation.

You can do this by using `-m`/`--modules` to explicitly define which modules
you want to run. Alternatively, use `-e`/`--exclude` to run all modules
except those listed.