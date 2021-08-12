# Configuring MultiQC

Whilst most MultiQC settings can be specified on the command line, MultiQC
is also able to parse system-wide and personal config files. At run time,
it collects the configuration settings from the following places in this order
(overwriting at each step if a conflicting config variable is found):

1. Hardcoded defaults in MultiQC code
2. System-wide config in `<installation_dir>/multiqc_config.yaml`
   - Manual installations only, not `pip` or `conda`
3. User config in `~/.multiqc_config.yaml`
4. File path set in environment variable `MULTIQC_CONFIG_PATH`
   - For example, define this in your `~/.bashrc` file and keep the file anywhere you like
5. Config file in the current working directory: `multiqc_config.yaml`
6. Config file paths specified in the command with `--config` / `-c`
   - You can specify multiple files like this, they can have any filename.
7. Command line config (`--cl_config`)
8. Specific command line options (_e.g._ `--force`)

## Sample name cleaning

MultiQC typically generates sample names by taking the input or log file name,
and 'cleaning' it. To do this, it uses the `fn_clean_exts` settings and looks
for any matches. If it finds any matches, everything to the right is removed.
For example, consider the following config:

```yaml
fn_clean_exts:
  - ".gz"
  - ".fastq"
```

This would make the following sample names:

```txt
mysample.fastq.gz  ->  mysample
secondsample.fastq.gz_trimming_log.txt  ->  secondsample
thirdsample.fastq_aligned.sam.gz  ->  thirdsample
```

There is also a config list called `fn_clean_trim` which just removes
strings if they are present at the start or end of the sample name.

Usually you don't want to overwrite the defaults (though you can).
Instead, add to the special variable names `extra_fn_clean_exts`
and `extra_fn_clean_trim`:

```yaml
extra_fn_clean_exts:
  - ".myformat"
  - "_processedFile"
extra_fn_clean_trim:
  - "#"
  - ".myext"
```

### Other search types

File name cleaning can also take strings to remove (instead of removing with truncation).
Also regex strings can be supplied to match patterns and remove or keep matching substrings.

#### `truncate` (default)

If you just supply a string, the default behavior is similar to "trim". The filename will be truncated beginning with the matching string.

```yaml
extra_fn_clean_exts:
  - ".fastq"
```

The above is equivalent to the more explicit:

```yaml
extra_fn_clean_exts:
  - type: "truncate"
    pattern: ".fastq"
```

This rule would produce the following sample names:

```txt
mysample.fastq.gz  ->  mysample
thirdsample.fastq_aligned.sam.gz  ->  thirdsample
```

#### `remove` (formerly `replace`)

The `remove` type allows you to remove the exact match from the filename.

```yaml
extra_fn_clean_exts:
  - type: remove
    pattern: .sorted
```

This rule would produce the following sample names:

```txt
secondsample.sorted.deduplicated  ->  secondsample.deduplicated
```

#### `regex`

You can also remove a substring with a regular expression. Here's a [good resource](https://regex101.com/) to interactively try it out.

```yaml
extra_fn_clean_exts:
  - type: regex
    pattern: "^processed."
```

This rule would produce the following sample names:

```txt
processed.thirdsample.processed  ->  thirdsample.processed
```

#### `regex_keep`

If you'd rather like to _keep_ the match of a regular expression you can use the `regex_keep` type. This simplifies things if you can e.g. directly target samples names.

```yaml
extra_fn_clean_exts:
  - type: regex_keep
    pattern: "[A-Z]{3}[1-9]{2}"
```

This rule would produce the following sample names:

```txt
merged.recalibrated.XZY97.alignment.bam  ->  XZY97
```

#### `module`

This key will tell MultiQC to only apply the pattern to a specific MultiQC module.
This should be a string that matches the module's `anchor` - the `#module` bit when you click the main module heading in the sidebar (remove the `#`).

For example, to truncate all sample names to 5 characters for just Kallisto:

```yaml
extra_fn_clean_exts:
  - type: regex_keep
    pattern: "^.{5}"
    module: kallisto
```

You can also supply a list of multiple module anchors if you wish:

```yaml
extra_fn_clean_exts:
  - type: regex_keep
    pattern: "^.{5}"
    module:
      - kallisto
      - cutadapt
```

### Clashing sample names

This process of cleaning sample names can sometimes result in exact duplicates.
A duplicate sample name will overwrite previous results. Warnings showing these events
can be seen with verbose logging using the `--verbose`/`-v` flag, or in `multiqc_data/multiqc.log`.

Problems caused by this will typically be discovered be fewer results than expected. If you're ever
unsure about where the data from results within MultiQC reports come from, have a look at
`multiqc_data/multiqc_sources.txt`, which lists the path to the file used for every section
of the report.

#### Directory names

One scenario where clashing names can occur is when the same file is processed in different directories.
For example, if `sample_1.fastq` is processed with four sets of parameters in four different
directories, they will all have the same name - `sample_1`. Only the last will be shown.
If the directories are different, this can be avoided with the `--dirs`/`-d` flag.

For example, given the following files:

```txt
├── analysis_1
│   └── sample_1.fastq.gz.aligned.log
├── analysis_2
│   └── sample_1.fastq.gz.aligned.log
└── analysis_3
    └── sample_1.fastq.gz.aligned.log
```

Running `multiqc -d .` will give the following sample names:

```txt
analysis_1 | sample_1
analysis_2 | sample_1
analysis_3 | sample_1
```

#### Filename truncation

If the problem is with filename truncation, you can also use the `--fullnames`/`-s` flag,
which disables all sample name cleaning. For example:

```txt
├── sample_1.fastq.gz.aligned.log
└── sample_1.fastq.gz.subsampled.fastq.gz.aligned.log
```

Running `multiqc -s .` will give the following sample names:

```txt
sample_1.fastq.gz.aligned.log
sample_1.fastq.gz.subsampled.fastq.gz.aligned.log
```

You can turn off sample name cleaning permanently by setting
`fn_clean_sample_names` to `false` in your config file.

## Module search patterns

Many bioinformatics tools have standard output formats, filenames and other
signatures. MultiQC uses these to find output; for example, the FastQC module
looks for files that end in `_fastqc.zip`.

This works well most of the time, until someone has an automated processing
pipeline that renames things. For this reason, as of version v0.3.2 of MultiQC,
the file search patterns are loaded as part of the main config. This means that
they can be overwritten in `<installation_dir>/multiqc_config.yaml` or
`~/.multiqc_config.yaml`. So if you always rename your `_fastqc.zip` files to
`_qccheck.zip`, MultiQC can still work.

To see the default search patterns, check a given module in the MultiQC documentation.
Each module has its search patterns listed beneath any free-text docs.
Alternatively, see the [`search_patterns.yaml`](https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml)
file in the MultiQC source code. Copy the section for the program that you want to modify and paste this
into your config file. Make sure you make it part of a dictionary called `sp`
as follows:

```yaml
sp:
  mqc_module:
    fn: _mysearch.txt
```

Search patterns can specify a filename match (`fn`) or a file contents
match (`contents`), as well as a number of additional search keys.
See [below](#step-1-find-log-files) for the full reference.

## Using log filenames as sample names

A substantial number of MultiQC modules take the sample name identifiers that you
see in the report from the file contents - typically the filename of the input file.
This is because log files can often be called things like `mytool.log` or even concatenated.
Using the input filename used by the tool is typically safer and more consistent across modules.

However, sometimes this does not work well. For example, if the input filename is not
relevant (eg. using a temporary file or FIFO, process substitution or stdin etc.).
In these cases your log files may have useful filenames but MultiQC will not be using them.

To force MultiQC to use the log filename as the sample identifier, you can use the
`--fn_as_s_name` command line flag or set the `use_filename_as_sample_name`:

```yaml
use_filename_as_sample_name: true
```

This affects all modules and all search patterns. If you want to limit this to just
one or more specific search patterns, you can do by giving a list:

```yaml
use_filename_as_sample_name:
  - cutadapt
  - picard/gcbias
  - picard/markdups
```

Note that this should be the search pattern key (see above) and not just the module name.
This is because some modules search for multiple files.

The log filename will still be cleaned. To use the raw log filenames,
combine with the `--fullnames`/`-s` flag or `fn_clean_sample_names` config option described above.

## Ignoring Files

MultiQC begins by indexing all of the files that you specified and building a list
of the ones it will use. You can specify files and directories to skip on the command
line using `-x`/`--ignore`, or for more permanent memory, with the following config file
options: `fn_ignore_files`, `fn_ignore_dirs` and `fn_ignore_paths` (the command line
option simply adds to all of these).

For example, given the following files:

```txt
├── analysis_1
│   └── sample_1.fastq.gz.aligned.log
├── analysis_2
│   └── sample_1.fastq.gz.aligned.log
└── analysis_3
    └── sample_1.fastq.gz.aligned.log
```

You could specify the following relevant config options:

```yaml
fn_ignore_files:
  - "*.log"
fn_ignore_dirs:
  - "analysis_1"
  - "analysis_2"
fn_ignore_paths:
  - "*/analysis_*/sample_1*"
```

Note that the searched file paths will usually be relative to the working
directory and can be highly variable, so you'll typically want to start patterns
with a `*` to match any preceding directory structure.

## Ignoring samples

Some modules get sample names from the contents of the file and not the filename
(for example, `stdout` logs can contain multiple samples). You can skip samples
by their resolved sample names (after cleaning) with two config options:
`sample_names_ignore` and `sample_names_ignore_re`. The first takes a list of
strings to be used for glob pattern matching (same behaviour as the command line
option `--ignore-samples`), the latter takes a list of regex patterns. For example:

```yaml
sample_names_ignore:
  - "SRR*"
sample_names_ignore_re:
  - '^SR{2}\d{7}_1$'
```

## Large sample numbers

MultiQC has been written with the intention of being used for any number of samples.
This means that it _should_ work well with 6 samples or 6000. Very large sample numbers
are becoming increasingly common, for example with single cell data.

Producing reports with data from many hundreds or thousands of samples provides some
challenges, both technically and also in terms of data visualisation and report usability.

### Disabling on-load plotting

One problem with large reports is that the browser can hang when the report is first loaded.
This is because it loading and processing the data for all plots at once. To mitigate this,
large reports may show plots as grey boxes with a _"Show Plot"_ button. Clicking this will
render the plot as normal and prevents the browser from trying to do everything at once.

By default this behaviour kicks in when a plot has 50 samples or more. This can be customised
by changing the `num_datasets_plot_limit` config option.

### Flat / interactive plots

Reports with many samples start to need a lot of data for plots. This results in inconvenient
report file sizes (can be 100s of megabytes) and worse, web browser crashes. To allow MultiQC
to scale to these sample numbers, most plot types have two plotting functions in the code base -
interactive (using HighCharts) and flat (rendered with MatPlotLib). Flat plots take up the
same disk space irrespective of sample number and do not consume excessive resources to display.

By default, MultiQC generates flat plots when there are 100 or more samples. This cutoff
can be changed by changing the `plots_flat_numseries` config option. This behaviour can also
be changed by running MultiQC with the `--flat` / `--interactive` command line options or by
setting the `plots_force_flat` / `plots_force_interactive` config options to `True`.

### Tables / Beeswarm plots

Report tables with thousands of samples (table rows) can quickly become impossible to use.
To avoid this, tables with large numbers of rows are instead plotted as a Beeswarm plot
(aka. a strip chart / jitter plot). These plots have fixed dimensions with any number
of samples. Hovering on a dot will highlight the same sample in other rows.

By default, MultiQC starts using beeswarm plots when a table has 500 rows or more. This
can be changed by setting the `max_table_rows` config option.

## Coloured log output

As of MultiQC version 1.8, log output is coloured using the [coloredlogs](https://pypi.org/project/coloredlogs/)
Python package. The code attempts to detect if the logs on the terminal are being redirected to a file
or piped to another tool and will disable colours if so. If the colours annoy you or you're ending
up with weird characters in your MultiQC output, you can disable this feature with the command line
flag `--no-ansi`. Sadly it's not possible to set this in a config file, as the logger is initilised
before configs are loaded.

## Command-line config

Sometimes it's useful to specify a single small config option just once, where creating
a config file for the occasion may be overkill. In these cases you can use the
`--cl_config` option to supply additional config values on the command line.

Config variables should be given as a YAML string. You will usually need to enclose
this in quotes. If MultiQC is unable to understand your config you will get an error message
saying `Could not parse command line config`.

As an example, the following command configures the coverage levels to use for the
Qualimap module: _(as [described in the docs](http://multiqc.info/docs/#qualimap))_

```bash
multiqc ./datadir --cl_config "qualimap_config: { general_stats_coverage: [20,40,200] }"
```

## Optimising run-time

Usually, MultiQC run time is fairly insignificant - in the order of seconds.
Unless you are running MultiQC on many thousands of analysis files, the optimisations
described below will have limited practical benefit.

In other words, if you're running with 15 RNAseq samples, you may as well save yourself
some time and stick with the defaults.

### Profile your MultiQC run time

As of version 1.9, MultiQC has a command line option to profile what it spends its time
doing: `--profile-runtime` (`config.profile_runtime`). Whilst you're working with writing
your pipeline / setting up your analysis, you can specify and MultiQC will add a section
to the bottom of your report describing how much time it spent searching files and
what it did with those files. You'll also get a breakdown in the command-line log
of how long the different steps of MultiQC execution took:

```txt
[INFO   ]         multiqc : MultiQC complete
[INFO   ]         multiqc : Run took 35.28 seconds
[INFO   ]         multiqc :  - 31.01s: Searching files
[INFO   ]         multiqc :  - 1.75s: Running modules
[INFO   ]         multiqc :  - 0.96s: Compressing report data
[INFO   ]         multiqc : For more information, see the 'Run Time' section in multiqc_report.html
```

If MultiQC is finishing in a few seconds or minutes, you probably don't need to do anything.
If you are working with huge numbers of files then it may be worth looking into these
results to see if you can speed up MultiQC. The documentation below explains how to do this.

### Be picky with which modules are run

Probably the easiest way to speed up MultiQC is to only use the modules that you
know you have files for. MultiQC supports a _lot_ of different tools and searches
for matching files for all of them every time you run it.

You can do this with the `-m` / `--module` flag (can be repeated) or in a MultiQC
config file by using `config.module_order`. See [Order of modules](#order-of-modules).

### Optimise file search patterns

Secondly, think about customising the search patterns of the slowest searches.

As an example, logs from Picard are published to `STDOUT` and so can have any file name.
Some people concatenate logs, so the contents can be anywhere in the file and the files
must also be searched by subsequent tools in case they contain multiple outputs.
If you know that all of your Picard MarkDuplicate log files have the filename
`mysamplename_markduplicates.log` then you can safely customise that search pattern
with the following MultiQC config:

```yaml
sp:
  picard/markdups:
    fn: "*_markduplicates.log"
```

If you know that this is the only type of Picard output that you're interested in,
you can also change all of the other Picard search patterns to use `skip: True`:

```yaml
sp:
  picard/markdups:
    fn: "*_markduplicates.log"
  picard/alignment_metrics:
    skip: true
  picard/basedistributionbycycle:
    skip: true
  picard/gcbias:
    skip: true
  picard/hsmetrics:
    skip: true
  picard/insertsize:
    skip: true
  picard/oxogmetrics:
    skip: true
  picard/pcr_metrics:
    skip: true
  picard/quality_by_cycle:
    skip: true
  picard/quality_score_distribution:
    skip: true
  picard/quality_yield_metrics:
    skip: true
  picard/rnaseqmetrics:
    skip: true
  picard/rrbs_metrics:
    skip: true
  picard/sam_file_validation:
    skip: true
  picard/variant_calling_metrics:
    skip: true
  picard/wgs_metrics:
    skip: true
```

This can speed up execution a bit if you really want to squeeze that running time.
The [MultiQC Modules documentation](#multiqc-modules) shows the search patterns for every module.

> Note that it's only worth using `skip: true` on search patterns if you want to use one from a module that has several.
> Usually it's better to just [specify which modules you want to run](#be-picky-with-which-modules-are-run) instead.

### Force interactive plots

One step that can take some time is running MatPlotLib to generate static-image plots
(see [Flat / interactive plots](#flat--interactive-plots)).
You can force MultiQC to skip this and only use interactive plots by using the `--interactive`
command line option (`config.plots_force_interactive`).

This approach is **not recommended if you have a very large number of samples**, as this can
produce a huge report file with all of the embedded plot data and crash your browser when opening it.
If you are running MultiQC for the `multiqc_data` folder and never intend to look at the report, it
speed things up though.
