---
title: Configuration
description: Settings to tweak how MultiQC works
---

# Configuring MultiQC

Whilst most MultiQC settings can be specified on the command line, MultiQC
is also able to parse system-wide and personal config files. At run time,
it collects the configuration settings from the following places in this order
(overwriting at each step if a conflicting config variable is found):

1. Hardcoded defaults in MultiQC code
1. System-wide config in `<installation_dir>/multiqc_config.yaml`
   - Manual installations only, not `pip` or `conda`
1. User config in `$XDG_CONFIG_HOME/multiqc_config.yaml` (or `~/.config/multiqc_config.yaml` if `$XDG_CONFIG_HOME` is not set)
1. User config in `~/.multiqc_config.yaml`
1. File path set in environment variable `$MULTIQC_CONFIG_PATH`
   - For example, define this in your `~/.bashrc` file and keep the file anywhere you like
1. Environment variables prefixed with `MULTIQC_`
   - For example, `$MULTIQC_TITLE`, `$MULTIQC_TEMPLATE` - see [docs below](#config-with-environment-variables)
1. Config file in the current working directory: `multiqc_config.yaml`
1. Config file paths specified in the command with `--config` / `-c`
   - You can specify multiple files like this, they can have any filename.
1. Command line config (`--cl-config`)
1. Specific command line options (_e.g._ `--force`)

## Sample name cleaning

MultiQC typically generates sample names by taking the input or log file name,
and 'cleaning' it.

### Cleaning extensions

To do this, it uses the `fn_clean_exts` settings and looks for any matches. If it finds any matches, everything to the right is removed.

:::info{title=Example}

```yaml
fn_clean_exts:
  - ".gz"
  - ".fastq"
```

| Input                                    | Cleaned sample name |
| ---------------------------------------- | ------------------- |
| `mysample.fastq.gz`                      | `mysample`          |
| `secondsample.fastq.gz_trimming_log.txt` | `secondsample`      |
| `thirdsample.vcf.gz_report.txt`          | `thirdsample.vcf`   |

:::

To add to the MultiQC defaults instead of overwriting them, use `extra_fn_clean_exts`:

```yaml
extra_fn_clean_exts:
  - ".myformat"
  - "_processedFile"
```

### Trimming extensions

To remove a substring only if it is at the start or end of a sample name, rather than trimming it and everything after it, use `fn_clean_trim`.

:::info{title=Example}

```yaml
fn_clean_trim:
  - ".fastq.gz"
  - "_report.txt"
```

| Input                                    | Cleaned sample name                      |
| ---------------------------------------- | ---------------------------------------- |
| `mysample.fastq.gz`                      | `mysample`                               |
| `secondsample.fastq.gz_trimming_log.txt` | `secondsample.fastq.gz_trimming_log.txt` |
| `thirdsample.vcf.gz_report.txt`          | `thirdsample.vcf.gz`                     |

:::

Again, to add to the MultiQC defaults instead of overwriting them, use `extra_fn_clean_trim`:

```yaml
extra_fn_clean_trim:
  - "#"
  - ".myext"
```

### Other search types

If needed, you can specify different string matching methods to `fn_clean_exts` and `extra_fn_clean_exts` for more complex sample name cleaning:

#### `truncate` (default)

This is the default method as described above. The two examples below are equivalent:

```yaml
extra_fn_clean_exts:
  - ".fastq"
```

```yaml
extra_fn_clean_exts:
  - type: "truncate"
    pattern: ".fastq"
```

#### `remove`

The `remove` type allows you to remove an exact match from the filename.
This includes removing a substring within the middle of a sample name.

:::info{title=Example}

```yaml
extra_fn_clean_exts:
  - type: remove
    pattern: .sorted
```

| Input                               | Cleaned sample name         |
| ----------------------------------- | --------------------------- |
|  `secondsample.sorted.deduplicated` | `secondsample.deduplicated` |

:::

#### `regex`

You can also remove a substring with a regular expression.
A useful website to work with writing regexes is [regex101.com](https://regex101.com).

:::info{title=Example}

```yaml
extra_fn_clean_exts:
  - type: regex
    pattern: "^processed."
```

| Input                             | Cleaned sample name     |
| --------------------------------- | ----------------------- |
| `processed.thirdsample.processed` | `thirdsample.processed` |

:::

#### `regex_keep`

If you'd rather like to only _keep_ the match of a regular expression and discard everything else, you can use the `regex_keep` type.
This is particularly useful if you have predictable sample names or identifiers.

:::info{title=Example}

```yaml
extra_fn_clean_exts:
  - type: regex_keep
    pattern: "[A-Z]{3}[1-9]{2}"
```

| Input                                     | Cleaned sample name |
| ----------------------------------------- | ------------------- |
| `merged.recalibrated.XZY97.alignment.bam` | `XZY97`             |

:::

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

```
├── analysis_1
│   └── sample_1.fastq.gz.aligned.log
├── analysis_2
│   └── sample_1.fastq.gz.aligned.log
└── analysis_3
    └── sample_1.fastq.gz.aligned.log
```

Running `multiqc -d .` will give the following sample names:

```
analysis_1 | sample_1
analysis_2 | sample_1
analysis_3 | sample_1
```

#### Filename truncation

If the problem is with filename truncation, you can also use the `--fullnames`/`-s` flag,
which disables all sample name cleaning. For example:

```
├── sample_1.fastq.gz.aligned.log
└── sample_1.fastq.gz.subsampled.fastq.gz.aligned.log
```

Running `multiqc -s .` will give the following sample names:

```
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
Alternatively, see the [`search_patterns.yaml`](https://github.com/MultiQC/MultiQC/blob/main/multiqc/search_patterns.yaml)
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
See [below](../development/modules.md#step-1---find-log-files) for the full reference.

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

Note that this should be the search pattern key and not just the module name.
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

```
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
by changing the `plots_num_samples_do_not_automatically_load` config option.

### Flat / interactive plots

Reports with many samples start to need a lot of data for plots. This results in inconvenient
report file sizes (can be 100s of megabytes) and worse, web browser crashes. To allow MultiQC
to scale to these sample numbers, most plot types have two plotting methods in the code base -
interactive and flat. Flat plots take up the same disk space irrespective of sample number
and do not consume excessive resources to display.

By default, MultiQC generates flat plots when there are 1000 or more samples. This cutoff
can be changed by changing the `plots_flat_numseries` config option. This behaviour can also
be changed by running MultiQC with the `--flat` / `--interactive` command line options or by
setting the `plots_force_flat` / `plots_force_interactive` config options to `True`.

### Tables / violin plots

Report tables with thousands of samples (table rows) can quickly become impossible to use.
To avoid this, tables with large numbers of rows are instead plotted as a violin plot.
These plots have fixed dimensions with any number of samples, and can be helpful
to see the data distribution of each table column. By default, MultiQC starts using
violin plots when a table has 500 rows or more. This can be changed by setting the
`max_table_rows` config option.

There are also interactive dots for separate samples that can be hovered to show
sample name and highlight this sample in other rows. For efficiency, if the number
of samples is above `violin_min_threshold_outliers` (default value 100), only dots
for outliers within the distribution are shown, and for more than
`violin_min_threshold_no_points` (1000) samples, only the violins without points are
rendered. When the number of samples is above `violin_downsample_after` (2000),
the underlying violin data itself is downsampled to keep the interactive reports efficient.

## Coloured log output

As of MultiQC version 1.8, log output is coloured using the [coloredlogs](https://pypi.org/project/coloredlogs/)
Python package. The code attempts to detect if the logs on the terminal are being redirected to a file
or piped to another tool and will disable colours if so. If the colours annoy you or you're ending
up with weird characters in your MultiQC output, you can disable this feature with the command line
flag `--no-ansi`. Sadly it's not possible to set this in a config file, as the logger is initilised
before configs are loaded.

## Checks for new versions

When MultiQC runs it automatically checks to see if there is a new version available to download.
If so, a log message is printed at the top of the run saying where to download it
(_MultiQC Version v0.6 now available!_).
This helps people stay up to date and reduces the number of bug reports that are
due to outdated MultiQC versions.

The timeout for the version check is set to 5 seconds, so if you're running offline it should
fail silently and add negligable run time.
However, if you prefer you can explicitly disable the version check by adding
`no_version_check: true` to your MultiQC config.

The check is done with the main [MultiQC API](https://api.multiqc.info/)
(see [source code](https://github.com/MultiQC/api.multiqc.info)).
The only statistics that are collected are the number of checks and a handful of metrics about
the running environment of MultiQC, such as the Python version and installation method
(see [source code](https://github.com/MultiQC/MultiQC/blob/06faefd772ade811c3f9968d8db6106bd14eb57a/multiqc/core/version_check.py#L24-L34)).
No identifiable information (such as IP address) is stored.

## Command-line config

Sometimes it's useful to specify a single small config option just once, where creating
a config file for the occasion may be overkill. In these cases you can use the
`--cl-config` option to supply additional config values on the command line.

Config variables should be given as a YAML string. You will usually need to enclose
this in quotes. If MultiQC is unable to understand your config you will get an error message
saying `Could not parse command line config`.

As an example, the following command configures the coverage levels to use for the
Qualimap module: _(as [described in the docs](https://multiqc.info/modules/qualimap/))_

```bash
multiqc ./datadir --cl-config "qualimap_config: { general_stats_coverage: [20,40,200] }"
```

## Config with environment variables

Config parameters can be set through environment variables prefixed with `MULTIQC_`.
For example, setting:

```bash
export MULTIQC_TITLE="My report"
export MULTIQC_FILESEARCH_LINES_LIMIT=10
```

Is equivalent to setting these in YAML:

```yaml
title: "My report"
filesearch_lines_limit: 10
```

:::tip
Some variables such as `title` can be also directly set through the command line
options: `--title "My report"`. For a list of all parameters, run `multiqc --help`.
:::

Note that it is _not_ possible to set nested config parameters through environment
variables, such as those that expect lists or dicts as values (e.g. `fn_clean_exts`).

## Referencing environment variables in YAML configs

It is also supported to interpolate environment variables with config files. For
example, if you have a config file `multiqc_config.yaml` with the following content:

```yaml
title: !ENV "${TITLE}"
report_header_info:
  - Contact E-mail: !ENV "${NAME:info}@${DOMAIN:example.com}"
```

And you have the following environment variables set:

```bash
export TITLE="My report"
export NAME="John"
```

You will get the following report header:

```
*My report*
Contact E-mail: John@example.com
```

See that the `$DOMAIN` environment variable was not set, and the default value
`"example.com"` is used instead.

For more details on environment variable interpolation, refer to the documentation
of [pyaml_env](https://github.com/mkaranasou/pyaml_env), that is used by MultiQC
internally to process user YAML files.

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

```
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
config file by using `config.module_order`. See [Order of modules](../reports/customisation.md#order-of-modules).

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
The [MultiQC Modules documentation](../development/modules.md) shows the search patterns for every module.

:::tip
Note that it's only worth using `skip: true` on search patterns if you want to use one from a module that has several.
Usually it's better to just [specify which modules you want to run](#be-picky-with-which-modules-are-run) instead.
:::

### Force interactive plots

One step that can take some time is generating static-image plots
(see [Flat / interactive plots](https://multiqc.info/docs/development/plots/#interactive--flat-image-plots)).
You can force MultiQC to skip this and only use interactive plots by using the `--interactive`
command line option (`config.plots_force_interactive`).

This approach is **not recommended if you have a very large number of samples**, as this can
produce a huge report file with all of the embedded plot data and crash your browser when opening it.
If you are running MultiQC for the `multiqc_data` folder and never intend to look at the report, it
speed things up though.

### Skip the report if you don't need it

If you're running MultiQC just to get parsed data / exported plots (`multiqc_data`) or the output for MegaQC
and don't actually need the report, you can skip it with `--no-report`.
This prevents any HTML report from being generated, including the data compression step that precedes it.
This can cut a few seconds off the MultiQC execution time.

## Custom CSS files

MultiQC generates HTML reports. You can include custom CSS in your final report if you wish.
Simply add CSS files to the `custom_css_files` config option:

```yaml
custom_css_files:
  - myfile.css
```

Or pass `--custom-css-file` (can be specified multiple times) and MultiQC will include
them in the final report HTML.
