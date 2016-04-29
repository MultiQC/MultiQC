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
section of the docs (`usage.md`) explains this in more detail.

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

## Error messages about mkl trial mode / licences
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
can find instructions [here](http://multiqc.info/docs/#installing-with-conda)
(`installation.md`).

## Locale Error Messages
A more obscure problem is an error from the MatPlotLib python library
about locale settings. Some strings (such as `en_SE`) aren't allowed.
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

You can fix this by changing the locale to something that will be
recognised. One way to do this is by adding these lines to your 
`.bashrc` in your home directory (or `.bash_profile`):
```bash
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```


