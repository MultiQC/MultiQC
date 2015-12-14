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
You can use a custom name for the report with the `-n`/`--name` parameter, or instruct
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

## Printing to stdout
If you would like to generate MultiQC reports on the fly, you can print the
output to standard out by specifying `-n stdout`. Note that the data directory
will not be generated and the template used must create stand-alone HTML reports.

## Parsed data directory
By default, MultiQC creates a directory alongside the report containing
tab-delimited files with the parsed data. This is useful for downstream
processing, especially if you're running MultiQC with very large numbers
of samples.

You can explicitly choose whether to produce this data by specifying either the
`--data-dir` or `--no-data-dir` command line flags. You can also set a default 
in your configuration settings (`make_data_dir`). Note that the data directory
is never produced when printing the MultiQC report to `stdout`.

## Choosing which modules to run
Sometimes, it's desirable to choose which MultiQC modules run. This could be
because you're only interested in one type of output and want to keep the
reports small. Or perhaps the output from one module is misleading in your
situation.

You can do this by using `-m`/`--modules` to explicitly define which modules
you want to run. Alternatively, use `-e`/`--exclude` to run all modules
except those listed.

## MultiQC config settings
Whilst most MultiQC settings can be specified on the command line, MultiQC
is also able to parse system-wide and personal config files. At run time,
it collects the configuration settings from the following places in this order
(overwriting at each step if a conflicting config variable is found):

1. Hardcoded defaults in MultiQC code
2. System-wide config in `<installation_dir>/multiqc_config.yaml`
3. User config in `~/.multiqc_config.yaml`
4. Command line options

You can find an example configuration file that comes with MultiQC, called
`multiqc_config.yaml.example`