---
title: Common Problems
description: Troubleshooting difficulties
---

# Troubleshooting

Hopefully MultiQC will be easy to use and run without any hitches. If you have
any problems, please do get in touch with the developer
([Phil Ewels](http://phil.ewels.co.uk)) by e-mail or by
[submitting an issue](https://github.com/MultiQC/MultiQC/issues/new) on github.
Before that, here are a few things previously encountered that may help...

## Not enough samples found

In this scenario, MultiQC finds _some_ logs for the bioinformatics tool
in question, but not all of your samples appear in the report. This is
the most common question I get regarding MultiQC operation.

Usually, this happens because sample names collide. This happens innocently
a lot - MultiQC overwrites previous results of the same name and you get
the last one seen in the report. You can see warnings about this by running
MultiQC in verbose mode with the `-v` flag, or looking at the generated log
file in `multiqc_data/multiqc.log`. If you are unsure about what log file
ended up in the report, look at `multiqc_data/multiqc_sources.txt` which
lists each source file used.

To solve this, try running MultiQC with the `-d` and `-s` flags.
The [Clashing sample names](../getting_started/config.md#clashing-sample-names)
section of the docs explains this in more detail.

### Big log files

Another reason that log files can be skipped is if the log filesize is very
large. For example, this could happen with very long concatenated standard out
files.

By default, MultiQC skips any file that is larger than 50MB to keep
execution fast. The verbose log output (`-v` or `multiqc_data/multiqc.log`) will
show you if files are being skipped with messages such as these:

```
[DEBUG  ]  Ignoring file as too large: filename.txt
```

You can configure the threshold and parse your files by changing the
`log_filesize_limit` config option. For example, to parse files up to 2GB in
size, add the following to your MultiQC config file:

```yaml
log_filesize_limit: 2000000000
```

### Long log files

When MultiQC runs, it first scans all supplied input files to create a shortlist
of matching files for each module. These filenames are then passed to that module
in the second phase of the run time to parse.

Files are matches using search patterns that check for either filenames or specific
strings within a given file. When doing the latter, MultiQC will only look in the first
1000 lines of each file by default.

This should be fine for nearly all cases, but if you're concatenating log files
or are otherwise unlucky, the required string may fall beyond this limit.
There is no log message for this, as the file isn't being skipped - it's just an
absence of detection.

You can configure the threshold used by changing `filesearch_lines_limit` in your MultiQC config.
For example, to load all of every file (as was the default prior to MultiQC v1.15)
just set it to a very high number:

```yaml
filesearch_lines_limit: 2000000000
```

This will slow down the initial file search but should otherwise be safe.

### Concatenated log files

In order to speed up the file search, once MultiQC matches a file to a specific
module it ignores it for all others. The only time this can really cause problems
is if you concenate log outputs from multiple analysis steps into a single file.

This "consumption" of files during the search can be disabled in the module search
pattern configuration, by setting `shared: true`.
It can also be set using the dedicated config option `filesearch_file_shared`, to avoid
needing to modify the search pattern.
The latter takes a list of module names, for example:

```yaml
filesearch_file_shared:
  - cutadapt
  - star
```

## No logs found for a tool

In this case, you have run a bioinformatics tool and have some log files in
a directory. When you run MultiQC with that directory, it finds nothing
for the tool in question.

There are a couple of things you can check here:

1. Is the tool definitely
   [supported by MultiQC](https://github.com/MultiQC/MultiQC)? If not, why
   not [open an issue](https://github.com/MultiQC/MultiQC/issues/new) to
   request it!
2. Did your bioinformatics tool definitely run properly? I've spent quite
   a bit of time debugging MultiQC modules only to realise that the output
   files from the tool were empty or incomplete. If your data is missing,
   take a look and the raw files and make sure that there's something to see!
3. Did you make sure that the logs you're trying to run MultiQC on are the ones
   expected by the module in question? Check the [module documentation](https://multiqc.info/modules/)
   and have a look at the [example data](https://github.com/MultiQC/test-data/tree/main/data/modules) used for CI tests.

If everything looks fine, then MultiQC probably needs extending to support
your data. Tools have different versions, different parameters and different
output formats that can confuse the parsing code.
Please [open an issue](https://github.com/MultiQC/MultiQC/issues/new) with
your log files and we can get it fixed.

## Locale Error Messages

Click is know to throw this error when the locale isn't set:

```python
# ..truncated traceback.. #
File "click/_unicodefun.py", line 118, in _verify_python3_env 'for mitigation steps.' + extra)

RuntimeError: Click will abort further execution because Python 3 was configured to use ASCII
as encoding for the environment.  Consult http://click.pocoo.org/python3/for mitigation steps.
```

You can fix this by setting your system locale to something that will be recognised.
One way to do this is by adding these lines to your `.bashrc` in your home directory
(or `.zshrc`/`.bash_profile`/etc.):

```bash
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```

Other locale strings are also fine, as long as the variables are set and valid.

## No space left on device

If you're running with very large datasets or have an unusually small temporary file
storage location, you may run in to the following error:

```python
  File "/usr/local/lib/python3.8/site-packages/multiqc/utils/util_functions.py", line 71, in write_data_file
    print( jsonstr.encode('utf-8', 'ignore').decode('utf-8'), file=f)
OSError: [Errno 28] No space left on device
```

This happens because MultiQC writes all output files to a temporary directory before
moving them to their final location (depending on configuration, MultiQC can change its output filenames
based on the data that it parses).

To solve this problem you can manually set the temp folder to another folder that has more space.
This is best done with an environment variable which is understood by the base Python installation, `TMPDIR`.

For example:

```bash
export TMPDIR=/path/to/my/tmp
```

Note that you will have to do this in each session where you run MultiQC. To automatically apply
it at every login, add the above line to your `~/.bashrc` file.

Be aware that setting this could have unforeseen consequences as it could affect the behaviour of other tools.
