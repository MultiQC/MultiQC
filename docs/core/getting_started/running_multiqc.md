---
title: Running MultiQC
description: Generating MultiQC reports from your data
---

# Running MultiQC

Once installed, just go to your analysis directory and run `multiqc`, followed
by a list of directories to search. At it's simplest, this can just be `.`
(the current working directory):

```bash
multiqc .
```

That's it! MultiQC will scan the specified directories and produce a report
based on details found in any log files that it recognises.

See [Using MultiQC Reports](../reports/reports.md) for more information about how
to use the generated report.

For a description of all command line parameters, run `multiqc --help`.

:::info
Every command-line flag mentioned on this page has a corresponding configuration variable that can be set in a MultiQC config YAML file.
This may be preferable if using a lot of options, or running in a pipeline.
:::

## Choosing where to scan

You can supply MultiQC with as many directories or files as you like. Above,
we supply `.` - just the current directory, but all of these would work too:

```bash
multiqc data/
multiqc data/ ../proj_one/analysis/ /tmp/results
multiqc data/*_fastqc.zip
multiqc data/sample_1*
```

If the `--ignore-symlinks` flag is set, MultiQC will ignore symlinked directories and files.

### Ignoring files

You can also ignore files or directories using the `-x`/`--ignore` option.
This can be specified multiple times and accepts glob patterns (eg. using the `*` and `?` wildcards).

:::warning
Glob patterns should be enclosed in quotes to prevent them being expanded by bash.
:::

The argument can match filenames, directory names and entire paths. For example:

```bash
multiqc . --ignore "file"
multiqc . --ignore "fileA" --ignore "fileB"
multiqc . --ignore "_R?.zip"
multiqc . --ignore "run_two/*"
multiqc . --ignore "*/run_three/*/fastqc/*_R2.zip"
```

Some modules get sample names from the contents of the file and not the filename
(for example, `stdout` logs can contain multiple samples). In this case, you can
skip samples by name instead:

```bash
multiqc . --ignore-samples "sample_3*"
```

These strings are matched using glob logic (`*` and `?` are wildcards).

All of these settings can be saved in a MultiQC config file so that you don't have
to type them on the command line for every run.

### File of search paths

If you have a large list of specific files, you can supply a file containing a list of file paths, one per row.
MultiQC will only search the listed files.

```bash
multiqc --file-list my_file_list.txt
```

## Renaming reports

The report is called `multiqc_report.html` by default. Tab-delimited data files
are created in `multiqc_data/`, containing additional information.
You can use a custom name for the report with the `-n`/`--filename` parameter, or instruct
MultiQC to create them in a subdirectory using the `-o`/`--outdir` parameter.

Note that different MultiQC templates may have different defaults.

## Overwriting existing reports

It's quite common to repeatedly create new reports as new analysis results are generated. Instead of manually deleting old reports, you can just specify the `-f`/`--force` parameter and MultiQC will overwrite any conflicting report filenames.

## Choosing which modules to run

Sometimes, it's desirable to choose which MultiQC modules run. This could be because you're only interested in one type of output and want to keep the reports small. Or perhaps the output from one module is misleading in your situation.

You can do this by using `-m`/`--modules` to explicitly define which modules you want to run. Alternatively, use `-e`/`--exclude` to run all modules _except_ those listed.

If an explicitly requested module couldn't find any expected input files, MultiQC will
just continue with other modules. You can change this behaviour and make MultiQC
strict about missing input by setting the `--require-log` flag.
If set, MultiQC will exit with an error and exit code `1` if any of the modules specified with `-m` did not produce a section in the report.

## Directory prefixes in sample names

Sometimes, the same samples may be processed in different ways. If MultiQC
finds log files with the same sample name, the previous data will be overwritten
(this can be inspected by running MultiQC with `-v`/`--verbose`).

To avoid this, run MultiQC with the `-d`/`--dirs` parameter. This will prefix every
sample name with the directory path for that log file. As such, sample names should
now be unique, and not overwrite one-another.

By default, `--dirs` will prepend the entire path to each sample name. You can choose
which directories are added with the `-dd`/`--dirs-depth` parameter. Set to a positive
integer to use that many directories at the end of the path. A negative integer takes
directories from the start of the path.

For example, show the full relative file path in the sample name:

```
$ multiqc -d .
# analysis_1 | results | type | sample_1 | file.log
# analysis_2 | results | type | sample_2 | file.log
# analysis_3 | results | type | sample_3 | file.log
```

Prepend just the last directory name:

```
$ multiqc -d -dd 1 .
# sample_1 | file.log
# sample_2 | file.log
# sample_3 | file.log
```

Prepend the first directory name:

```
$ multiqc -d -dd -1 .
# analysis_1 | file.log
# analysis_2 | file.log
# analysis_3 | file.log
```

## Printing to stdout

If you would like to generate MultiQC reports on the fly, you can print the output to standard out by specifying `-n stdout`.
The data directory will not be generated and the template used must create stand-alone HTML reports.

## Using different templates

MultiQC is built around a templating system. You can produce reports with
different styling by using the `-t`/`--template` option. The available templates
are listed with `multiqc --help`.

If you're interested in creating your own custom template, see the
[writing new templates](../development/templates.md) section.

## Parsed data directory

By default, MultiQC creates a directory alongside the report containing
tab-delimited files with the parsed data. This is useful for downstream
processing, especially if you're running MultiQC with very large numbers
of samples.

Typically, these files are tab-delimited tables. However, you can get `JSON`
or `YAML` output for easier downstream parsing by specifying `-k`/`--data-format`
on the command line or `data_format` in your configuration file.

You can also choose whether to produce the data by specifying either the
`--data-dir` or `--no-data-dir` command line flags or the `make_data_dir`
variable in your configuration file. Note that the data directory
is never produced when printing the MultiQC report to `stdout`.

To zip the data directory, use the `-z`/`--zip-data-dir` flag.

## Exporting Plots

In addition to the HTML report, it's also possible to get MultiQC to save
plots as standalone files. You can do this with the `-p`/`--export` command
line flag. By default, plots will be saved in a directory called `multiqc_plots`
as `.png`, `.svg` and `.pdf` files. Raw data for the plots are also saved to files.

You can instruct MultiQC to always do this by setting the `export_plots` config
option to `true`, though note that this will add a few seconds on to execution time.
The `plots_dir_name` changes the default directory name for plots and the
`export_plot_formats` specifies what file formats should be created (must be
supported by Plotly).

Note that not all plot types are yet supported, so you may find some plots are
missing.

:::note
You can always save static image versions of plots from within MultiQC reports, using the [Export toolbox](../reports/reports#exporting-plots) in the side bar.
:::

## PDF Reports

Whilst HTML is definitely the format of choice for MultiQC reports due to
the interactive features that it can offer, PDF files are an integral part
of some people's workflows. To try to accommodate this, MultiQC has a
`--pdf` command line flag which will try to create a PDF report for you.

:::danger
PDF export support for MultiQC can be difficult to use and disables many core MultiQC features and even some plots. It should only be used as a last resort.
:::

To generate PDFs, MultiQC uses the `simple` template. This uses flat plots,
has no navigation or toolbar and strips out all JavaScript. The resulting
HTML report is pretty basic, but this simplicity is helpful when generating
PDFs.

Once the report is generated MultiQC attempts to call [Pandoc](http://pandoc.org/),
a command line tool able to convert documents between different file formats.
**You must have Pandoc already installed for this to work**. If you don't have
Pandoc installed, you will get an error message that looks like this:

```
Error creating PDF - pandoc not found. Is it installed? http://pandoc.org/
```

Please note that Pandoc is a complex tool and has a number of its own dependencies
for PDF generation. Notably, it uses LaTeX / XeLaTeX which you must also have installed.
Please make sure that you have the latest version of Pandoc and
that it can successfully convert basic HTML files to PDF before reporting
and errors.

Error messages from Pandoc are piped through to the MultiQC log,
for example if the xelatex dependency is not installed you will see the following:

```
xelatex not found. Please select a different --pdf-engine or install xelatex
```

Note that not all plots have flat image equivalents, so
some will be missing (at time of writing: FastQC sequence content plot,
beeswarm dot plots, heatmaps).
