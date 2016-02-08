# Configuring MultiQC
Whilst most MultiQC settings can be specified on the command line, MultiQC
is also able to parse system-wide and personal config files. At run time,
it collects the configuration settings from the following places in this order
(overwriting at each step if a conflicting config variable is found):

1. Hardcoded defaults in MultiQC code
2. System-wide config in `<installation_dir>/multiqc_config.yaml`
3. User config in `~/.multiqc_config.yaml`
4. Command line options

You can find an example configuration file that comes with MultiQC, called
`multiqc_config.yaml.example` - hopefully this should be self explanatory
by way of the comments. Common changes are discussed in more detail below.

## Sample name cleaning
MultiQC typically generates sample names by taking the input or log file name,
and 'cleaning' it. To do this, it uses the `fn_clean_exts` settings and looks
for any matches. If it finds any matches, everything to the right is removed.
For example, consider the following config:
```yaml
fn_clean_exts:
    - .gz
    - .fastq
```
This would make the following sample names:
```
mysample.fastq.gz  ->  mysample
secondsample.fastq.gz_trimming_log.txt  ->  secondsample
thirdsample.fastq_aligned.sam.gz  ->  thirdsample
```
Note that this process of cleaning sample names can result in duplicate
names. A duplicate sample name will overwrite previous results. This can
be seen with verbose logging using the `-v` flag.


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

