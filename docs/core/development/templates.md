---
title: Writing new templates
description: Making MultiQC reports your own
---

# Writing New Templates

MultiQC is built around a templating system that uses the
[Jinja](http://jinja.pocoo.org/) python package. This makes it very
easy to create new report templates that fit your needs.

## Core or plugin

If your template could be of use to others, it would be great if you
could add it to the main MultiQC package. You can do this by creating a
fork of the [MultiQC GitHub repository](https://github.com/MultiQC/MultiQC),
adding your template and then creating a pull request to merge your changes
back to the main repository.

If it's very specific template, you can create a new Python package which
acts as a plugin. For more information about this, see the
[plugins documentation](plugins.md).

## Creating a template skeleton

For a new template to be recognised by MultiQC, it must be a python submodule
directory with a `__init__.py` file. This must be referenced in the `setup.py`
installation script as an
[entry point](http://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins).

You can see the bundled templates defined in this way:

```python
entry_points = {
    'multiqc.templates.v1': [
        'default = multiqc.templates.default',
        'simple = multiqc.templates.simple',
        'geo = multiqc.templates.geo',
    ]
}
```

Note that these entry points can point to any Python modules, so if you're
writing a plugin module you can specify your module name instead. Just make
sure that `multiqc.templates.v1` is the same.

Once you've added the entry point, remember to install the package again:

```bash
pip install -e .
```

Using `-e` tells `pip` to softlink the plugin files instead of
copying, so changes made whilst editing files will be reflected when you
run MultiQC.

The `__init__.py` files must define two variables - the path to the template
directory and the main jinja template file:

```python
template_dir = os.path.dirname(__file__)
base_fn = 'base.html'
```

## Child templates

The default MultiQC template contains a _lot_ of code. Importantly, it includes
1448 lines of custom JavaScript (at time of writing) which powers the plotting
and dynamic functions in the report. You probably don't want to rewrite all of
this for your template, so to make your life easier you can create a
_child template_.

To do this, add an extra variable to your template's `__init__.py`:

```python
template_parent = 'default'
```

This tells MultiQC to use the template files from the `default` template unless
a file with the same name is found in your child template. For instance, if
you just want to add your own logo in the header of the reports, you can create
your own `header.html` which will overwrite the default header.

Files within the default template have comments at the top explaining what
part of the report they generate.

## Extra init variables

There are a few extra variables that can be added to the `__init__.py` file
to change how the report is generated.

Setting `output_dir` instructs MultiQC to put the report and it's contents
into a subdirectory. Set the string to your desired name. Note that this will
be prefixed if `-p`/`--prefix` is set at run time.

Secondly, you can copy additional files with your report when it is generated.
This is usually used to copy required images or scripts with the report. These
should be a list of file or directory paths, relative to the `__init__.py` file.
Directory contents will be copied recursively.

You can also override config options in the template. For example, setting
the value of `config.plots_force_flat` can force the report to only have
static image plots.

```python
from multiqc.utils import config

output_subdir = 'multiqc_report'
copy_files = ['assets']
config.plots_force_flat = True
```

## Jinja template variables

There are a number of variables that you can use within your Jinja template.
Two namespaces are available - `report` and `config`. You can print these
using the Jinja curly brace syntax, _eg._ `{{ config.version }}`. See the
[Jinja2 documentation](http://jinja.pocoo.org/docs/dev/templates/) for more
information.

The default MultiQC template includes dependencies in the HTML so that the
report is standalone. If you would like to do the same, use the `include_file`
function. For example:

```html
<script>
  {
    {
      include_file("js/jquery.min.js");
    }
  }
</script>
<img src="data:image/png;base64,{{ include_file('img/logo.png', b64=True) }}" />
```

## Appendices

### Custom plotting functions

If you don't like the default plotting functions built into MultiQC, you
can write your own! If you create a callable variable in a template called
either `bargraph` or `linegraph`, MultiQC will use that instead. For example:

```python
def custom_linegraph(plotdata, pconfig):
    return '<h1>Awesome line graph here</h1>'
linegraph = custom_linegraph

def custom_bargraph(plotdata, plotseries, pconfig):
    return '<h1>Awesome bar graph here</h1>'
bargraph = custom_bargraph
```

These particular examples don't do very much, but hopefully you get the idea.
Note that you have to set the variable `linegraph` or `bargraph` to your function.
