# Configuring MultiQC
Whilst most MultiQC settings can be specified on the command line, MultiQC
is also able to parse system-wide and personal config files. At run time,
it collects the configuration settings from the following places in this order
(overwriting at each step if a conflicting config variable is found):

1. Hardcoded defaults in MultiQC code
2. System-wide config in `<installation_dir>/multiqc_config.yaml`
3. User config in `~/.multiqc_config.yaml`
4. Config file in the current working directory: `multiqc_config.yaml`
5. Command line options

You can find an example configuration file bundled with MultiQC, called
`multiqc_config.example.yaml` - hopefully this should be self explanatory
through the included comments. Common changes are discussed in more detail below.

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

Usually you don't want to overwrite the defaults (though you can).
Instead, add to the special variable name `extra_fn_clean_exts`:

```yaml
extra_fn_clean_exts:
    - '.myformat'
    - '_processedFile'
```

### Other search types
As of MultiQC v0.6, file name cleaning can also take strings to remove (instead of 
removing with truncation). Also regex strings can be supplied to match patterns and remove strings.

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

