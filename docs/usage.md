# Running MultiQC
Once installed, just go to your analysis directory and run `multiqc`, followed
by a list of directories to search. At it's simplest, this can just be `.`
(the current working directory):
```
multiqc .
```

That's it! MultiQC will scan the specified directories and produce a report
based on details found in any log files that it recognises.

See [Using MultiQC Reports](http://multiqc.info/docs/#using-multiqc-reports) for more information about how
to use the generated report.

For a description of all command line parameters, run `multiqc --help`.

## Choosing where to scan
You can supply MultiQC with as many directories or files as you like. Above,
we supply `.` - just the current directory, but all of these would work too:
```
multiqc data/
multiqc data/ ../proj_one/analysis/ /tmp/results
multiqc data/*_fastqc.zip
multiqc data/sample_1*
```

You can also ignore files using the `-x`/`--ignore` flag (can be specified multiple
times). This takes a string which it matches using glob expansion to filenames,
directory names and entire paths:
```
multiqc . --ignore *_R2*
multiqc . --ignore run_two/
multiqc . --ignore */run_three/*/fastqc/*_R2.zip
```

Some modules get sample names from the contents of the file and not the filename
(for example, `stdout` logs can contain multiple samples). In this case, you can
skip samples by name instead:
```
multiqc . --ignore-samples sample_3*
```
These strings are matched using glob logic (`*` and `?` are wildcards).

All of these settings can be saved in a MultiQC config file so that you don't have
to type them on the command line for every run.

Finally, you can supply a file containing a list of file paths, one per row.
MultiQC only search the listed files.
```
multiqc --file-list my_file_list.txt
```

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

By default, `--dirs` will prepend the entire path to each sample name. You can choose
which directories are added with the `-dd`/`--dirs-depth` parameter. Set to a positive
integer to use that many directories at the end of the path. A negative integer takes
directories from the start of the path.

For example:
```
$ multiqc -d .
# analysis_1 | results | type | sample_1 | file.log
# analysis_2 | results | type | sample_2 | file.log
# analysis_3 | results | type | sample_3 | file.log

$ multiqc -d -dd 1 .
# sample_1 | file.log
# sample_2 | file.log
# sample_3 | file.log

$ multiqc -d -dd -1 .
# analysis_1 | file.log
# analysis_2 | file.log
# analysis_3 | file.log
```

## Using different templates
MultiQC is built around a templating system. You can produce reports with
different styling by using the `-t`/`--template` option. The available templates
are listed with `multiqc --help`.

If you're interested in creating your own custom template, see the
[writing new templates](http://multiqc.info/docs/#writing-new-templates) section.

## PDF Reports
Whilst HTML is definitely the format of choice for MultiQC reports due to
the interactive features that it can offer, PDF files are an integral part
of some people's workflows. To try to accommodate this, MultiQC has a
`--pdf` command line flag which will try to create a PDF report for you.

To do this, MultiQC uses the `simple` template. This uses flat plots,
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

Please note that Pandoc is a complex tool and uses LaTeX / XeLaTeX for PDF
generation. Please make sure that you have the latest version of Pandoc and
that it can successfully convert basic HTML files to PDF before reporting
and errors. Also note that not all plots have flat image equivalents, so
some will be missing (at time of writing: FastQC sequence content plot,
beeswarm dot plots, heatmaps).

## Printing to stdout
If you would like to generate MultiQC reports on the fly, you can print the
output to standard out by specifying `-n stdout`. Note that the data directory
will not be generated and the template used must create stand-alone HTML reports.

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
plots as stand alone files. You can do this with the `-p`/`--export` command
line flag. By default, plots will be saved in a directory called `multiqc_plots`
as `.png`, `.svg` and `.pdf` files.

You can instruct MultiQC to always do this by setting the `export_plots` config
option to `true`, though note that this will add a few seconds on to execution time.
The `plots_dir_name` changes the default directory name for plots and the
`export_plot_formats` specifies what file formats should be created (must be
supported by MatPlotLib).

Note that not all plot types are yet supported, so you may find some plots are
missing.

> Note: You can always save static image versions of plots from within
> MultiQC reports, using the [Export toolbox](http://multiqc.info/docs/#export) in the side bar.

## Choosing which modules to run
Sometimes, it's desirable to choose which MultiQC modules run. This could be
because you're only interested in one type of output and want to keep the
reports small. Or perhaps the output from one module is misleading in your
situation.

You can do this by using `-m`/`--modules` to explicitly define which modules
you want to run. Alternatively, use `-e`/`--exclude` to run all modules
except those listed.
