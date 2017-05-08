# Configuring MultiQC
Whilst most MultiQC settings can be specified on the command line, MultiQC
is also able to parse system-wide and personal config files. At run time,
it collects the configuration settings from the following places in this order
(overwriting at each step if a conflicting config variable is found):

1. Hardcoded defaults in MultiQC code
2. System-wide config in `<installation_dir>/multiqc_config.yaml`
  * Manual installations only, not `pip` or `conda`
3. User config in `~/.multiqc_config.yaml`
4. File path set in environment variable `MULTIQC_CONFIG_PATH`
  * For example, define this in your `~/.bashrc` file and keep the file anywhere you like
5. Config file in the current working directory: `multiqc_config.yaml`
6. Command line options

You can find an example configuration file with the MultiQC source code, called
[`multiqc_config.example.yaml`](https://github.com/ewels/MultiQC/blob/master/multiqc_config_example.yaml).
If you installed MultiQC with `pip` or `conda` you won't have this file locally,
but you can find it on GitHub:
[github.com/ewels/MultiQC](https://github.com/ewels/MultiQC/blob/master/multiqc_config_example.yaml).

## Sample name cleaning
MultiQC typically generates sample names by taking the input or log file name,
and 'cleaning' it. To do this, it uses the `fn_clean_exts` settings and looks
for any matches. If it finds any matches, everything to the right is removed.
For example, consider the following config:
```yaml
fn_clean_exts:
    - '.gz'
    - '.fastq'
```
This would make the following sample names:
```
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
    - '.myformat'
    - '_processedFile'
extra_fn_clean_trim:
    - '#'
    - '.myext'
```

### Other search types
File name cleaning can also take strings to remove (instead of removing with truncation).
Also regex strings can be supplied to match patterns and remove strings.

Consider the following:
```yaml
extra_fn_clean_exts:
    - '.fastq'
    - type: 'replace'
      pattern: '.sorted'
    - type: 'regex'
      pattern: '^processed.'
```
This would make the following sample names:
```
mysample.fastq.gz  ->  mysample
secondsample.sorted.deduplicated.fastq.gz_processed.txt  ->  secondsample.deduplicated
processed.thirdsample.fastq_aligned.sam.gz  ->  thirdsample
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

To see the default search patterns, see the
[`search_patterns.yaml`](https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml)
file. Copy the section for the program that you want to modify and paste this
into your config file. Make sure you make it part of a dictionary called `sp`
as follows:

```yaml
sp:
    mqc_module:
        fn: _mysearch.txt
```

Search patterns can specify a filename match (`fn`) or a file contents
match (`contents`).

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
    - '*.log'
fn_ignore_dirs:
    - 'analysis_1'
    - 'analysis_2'
fn_ignore_paths:
    - '*/analysis_*/sample_1*'
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
    - 'SRR*'
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

