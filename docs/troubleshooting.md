# Troubleshooting

Hopefully MultiQC will be easy to use and run without any hitches. If you have
any problems, please do get in touch with the developer
([Phil Ewels](http://phil.ewels.co.uk)) by e-mail or by
[submitting an issue](https://github.com/ewels/MultiQC/issues/new) on github.
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
The [Clashing sample names](http://multiqc.info/docs/#clashing-sample-names)
section of the docs explains this in more detail.

### Big log files

Another reason that log files can be skipped is if the log filesize is very
large. For example, this could happen with very long concatenated standard out
files. By default, MultiQC skips any file that is larger than 10MB to keep
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

## No logs found for a tool
In this case, you have run a bioinformatics tool and have some log files in
a directory. When you run MultiQC with that directory, it finds nothing
for the tool in question.

There are a couple of things you can check here:

1. Is the tool definitely
   [supported by MultiQC](https://github.com/ewels/MultiQC)? If not, why
   not [open an issue](https://github.com/ewels/MultiQC/issues/new) to
   request it!
2. Did your bioinformatics tool definitely run properly? I've spent quite
   a bit of time debugging MultiQC modules only to realise that the output
   files from the tool were empty or incomplete. If your data is missing,
   take a look and the raw files and make sure that there's something to see!

If everything looks fine, then MultiQC probably needs extending to support
your data. Tools have different versions, different parameters and different
output formats that can confuse the parsing code.
Please [open an issue](https://github.com/ewels/MultiQC/issues/new) with
your log files and we can get it fixed.

## Error messages about mkl trial mode / licences
In this case you run MultiQC and get something like this:

```
$ multiqc .

Vendor:  Continuum Analytics, Inc.
Package: mkl
Message: trial mode EXPIRED 2 days ago

    You cannot run mkl without a license any longer.
    A license can be purchased it at: http://continuum.io
    We are sorry for any inconveniences.

    SHUTTING DOWN PYTHON INTERPRETER
```


The `mkl` library provides optimisations for `numpy`, a requirement of
`MatPlotLib`. Recent versions of Conda have a bundled version which should
come with a licence and remove the warning. See
[this page](https://docs.continuum.io/mkl-optimizations/index#dismissing-mkl-trial-warnings)
for more info. If you already have Conda installed you can get the updated
version by running:
```
conda remove mkl-rt
conda install -f mkl
```

Another way around it is to uninstall `mkl`. It seems that `numpy` works
without it fine:
```
$ conda remove --features mkl
```
Problem solved! See more
[here](http://stackoverflow.com/questions/25204021/anaconda-running-python-cannot-run-mkl-without-a-license) and
[here](https://www.continuum.io/blog/developer-blog/anaconda-25-release-now-mkl-optimizations).

If you're not using Conda, try installing MultiQC with that instead. You
can find instructions [here](http://multiqc.info/docs/#installing-with-conda).

## Locale Error Messages
Two MultiQC dependencies have been known to throw errors due to problems
with the Python locale settings, or rather the lack of those settings.

MatPlotLib can complain that some strings (such as `en_SE`) aren't allowed.
Running MultiQC gives the following error:
```bash
$ multiqc --version
```
```python
# ..long traceback.. #
 File "/sw/comp/python/2.7.6_milou/lib/python2.7/locale.py", line 443, in _parse_localename
   raise ValueError, 'unknown locale: %s' % localename
ValueError: unknown locale: UTF-8
```

Click can have a similar problem if the locale isn't set when using
Python 3. That generates an error that looks like this:

```python
# ..truncated traceback.. #
File "click/_unicodefun.py", line 118, in _verify_python3_env 'for mitigation steps.' + extra)

RuntimeError: Click will abort further execution because Python 3 was configured to use ASCII
as encoding for the environment.  Consult http://click.pocoo.org/python3/for mitigation steps.
```

You can fix both of these problems by changing your system locale
to something that will be recognised. One way to do this is by adding
these lines to your `.bashrc` in your home directory (or `.bash_profile`):

```bash
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```

Other locale strings are also fine, as long as the variables are set and valid.

