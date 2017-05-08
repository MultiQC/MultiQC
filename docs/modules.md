# Writing New Modules

## Introduction
Writing a new module can at first seem a daunting task. However, MultiQC
has been written _(and refactored)_ to provide a lot of functionality
as common functions.

Provided that you are familiar with writing Python and you have a read
through the guide below, you should be on your way in no time!

If you have any problems, feel free to contact the author - details
here: [@ewels](https://github.com/ewels)

## Core modules / plugins
New modules can either be written as part of MultiQC or in a stand-alone
plugin. If your module is for a publicly available tool, please add it
to the main program and contribute your code back when complete via a
pull request.

If your module is for something _very_ niche, which no-one else can use,
you can write it as part of a custom plugin. The process is almost identical,
though it keeps the code bases separate. For more information about this,
see the docs about _MultiQC Plugins_ below.

## Initial setup
### Submodule
MultiQC modules are Python submodules - as such, they need their own
directory in `/multiqc/` with an `__init__.py` file. The directory should
share its name with the module. To follow common practice, the module
code usually then goes in a separate python file (also with the same name)
which is then imported by `__init__.py`:
```python
from __future__ import absolute_import
from .modname import MultiqcModule
```

### Entry points
Once your submodule files are in place, you need to tell MultiQC that they
are available as an analysis module. This is done within `setup.py` using
[entry points](https://pythonhosted.org/setuptools/setuptools.html#dynamic-discovery-of-services-and-plugins).
In `setup.py` you will see some code that looks like this:
```python
entry_points = {
    'multiqc.modules.v1': [
        'bismark = multiqc.modules.bismark:MultiqcModule',
        [...]
    ]
}
```
Copy one of the existing module lines and change it to use your module name.
The order is irrelevant, so stick to alphabetical if in doubt.
Once this is done, you will need to update your installation of MultiQC:
```
python setup.py develop
```

### MultiQC config
So that MultiQC knows what order modules should be run in, you need to add
your module to the core config file.

In `multiqc/utils/config_defaults.yaml` you should see a list variable called `module_order`.
This contains the name of modules in order of precedence. Add your module here
in an appropriate position.

### Documentation
Next up, you need to create a documentation file for your module. The reason
for this is twofold: firstly, docs are important to help people to use, debug
and extend MultiQC (you're reading this, aren't you?). Secondly,
having the file there with the appropriate YAML front matter will make the
module show up on the [MultiQC homepage](http://multiqc.info) so that everyone
knows it exists. This process is automated once the file is added to the core
repository.

This docs file should be placed in `docs/modules/<your_module_name>.md` and
should have the following structure:

```
---
Name: Tool Name
URL: http://www.amazing-bfx-tool.com
Description: >
    This amazing tool does some really cool stuff. You can describe it
    here and split onto multiple lines if you want. Not too long though!
---

Your documentation goes here. Feel free to use markdown and write whatever
you think would be helpful. Please avoid using heading levels 1 to 3.
```

Make a reference to this in the YAML frontmatter at the top of
`docs/README.md` - this allows the website to find the file to build
the documentation.

### Readme and Changelog
Last but not least, remember to add your new module to the main `README.md`
file and `CHANGELOG.md`, so that people know that it's there. Feel free to
add your name to the list of credits at the bottom of the readme.

### MultiqcModule Class
If you've copied one of the other entry point statements, it will have
ended in `:MultiqcModule` - this tells MultiQC to try to execute a class or
function called `MultiqcModule`.

To use the helper functions bundled with MultiQC, you should extend this
class from `multiqc.modules.base_module.BaseMultiqcModule`. This will give
you access to a number of functions on the `self` namespace. For example:
```python
from multiqc.modules.base_module import BaseMultiqcModule

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='My Module', anchor='mymod',
        href="http://www.awesome_bioinfo.com/my_module",
        info="is an example analysis module used for writing documentation.")
```

Ok, that should be it! The `__init__()` function will now be executed every
time MultiQC runs. Try adding a `print("Hello World!")` statement and see
if it appears in the MultiQC logs at the appropriate time...

Note that the `__init__` variables are used to create the header, URL link,
analysis module credits and description in the report.

### Logging
Last thing - MultiQC modules have a standardised way of producing output,
so you shouldn't really use `print()` statements for your `Hello World` ;)

Instead, use the `logger` module as follows:
```python
import logging
log = logging.getLogger(__name__)
# Initialise your class and so on
log.info('Hello World!')
```

Log messages can come in a range of formats:

* `log.debug`
  * Thes only show if MultiQC is run in `-v`/`--verbose` mode
* `log.info`
  * For more important status updates
* `log.warning`
  * Alert user about problems that don't halt execution
* `log.error` and `log.critical`
  * Not often used, these are for show-stopping problems


## Step 1 - Find log files
The first thing that your module will need to do is to find analysis log
files. You can do this by searching for a filename fragment, or a string
within the file. It's possible to search for both (a match on either
will return the file) and also to have multiple strings possible.

First, add your default patterns to:
```
MULTIQC_ROOT/multiqc/utils/search_patterns.yaml
```

Each search has a yaml key (typically the name of your module), with one
or more search criteria. The following keys can be used:

* `fn`
  * A glob filename pattern, used with the Python [`fnmatch`](https://docs.python.org/2/library/fnmatch.html) function
* `fn_re`
  * A regex filename pattern
* `contents`
  * A string to match within the file contents (checked line by line)
* `contents_re`
  * A regex to match within the file contents (checked line by line)
  * NB: Regex must match entire line (add `.*` to start and end of pattern to avoid this)
* `num_lines`
  * The number of lines to search through for the `contents` string. Default: all lines.
* `shared`
  * By default, once a file has been assigned to a module it is not searched again. Specify `shared: true` when your file can be shared between multiple tools (for example, part of a `stdout` stream).
* `max_filesize`
  * Files larger than the `log_filesize_limit` config key (default: 10MB) are skipped. If you know your files will be smaller than this and need to search by contents, you can specify this value (in bytes) to skip any files smaller than this limit.

Please try to use `num_lines` and `max_filesize` where possible as they will speed up
MultiQC execution time.

For example, two typical modules could specify search patterns as follows:

```yaml
mymod:
    fn: '_myprogram.txt'
myothermod:
    contents: 'This is myprogram v1.3'
```

Note that if you want to find multiple different log files for a single module,
the convention is to use the module name and a forward slash separator _(this changed
in the v1.0 release and the slashes help backwards-compatibility)_:
```yaml
fastqc/data:
    fn: 'fastqc_data.txt'
fastqc/zip:
    fn: '_fastqc.zip'
```

You can also supply a list of different patterns for a single log file type if needed.
If any of the patterns are matched, the file will be returned:
```yaml
mymod:
    - fn: 'mylog.txt'
    - fn: 'different_fn.out'
```

You can use _AND_ logic by specifying keys within a single list item. For example:
```yaml
mymod:
    fn: 'mylog.txt'
    contents: 'mystring'
myothermod:
    - fn: 'different_fn.out'
      contents: 'This is myprogram v1.3'
    - fn: 'another.txt'
      contents: 'What are these files anyway?'
```
Here, a file must have the filename `mylog.txt` _and_ contain the string `mystring`.

Remember that users can overwrite these defaults in their own config files.
This is helpful as people have weird and wonderful processing pipelines with
their own conventions.

Once your strings are added, you can find files in your module with the
base function `self.find_log_files()`, using the key you set in the YAML:
```python
self.find_log_files('mymod')
```

This function yields a dictionary with various information about each matching
file. The `f` key contains the contents of the matching file:
```python
# Find all files for mymod
for myfile in self.find_log_files('mymod'):
    print( myfile['f'] )       # File contents
    print( myfile['s_name'] )  # Sample name (from cleaned filename)
    print( myfile['fn'] )      # Filename
    print( myfile['root'] )    # Directory file was in
```

If `filehandles=True` is specified, the `f` key contains a file handle
instead:
```python
for f in self.find_log_files('mymod', filehandles=True):
    # f['f'] is now a filehandle instead of contents
    for l in f['f']:
        print( l )
```
This is good if the file is large, as Python doesn't read the entire
file into memory in one go.

## Step 2 - Parse data from the input files
What most MultiQC modules do once they have found matching analysis files
is to pass the matched file contents to another function, responsible
for parsing the data from the file. How this parsing is done will depend
on the format of the log file and the type of data being read. See below
for a basic example, based loosly on the preseq module:

```python
class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # [...]
        self.mod_data = dict()
        for f in self.find_log_files('mymod'):
            self.mod_data[f['s_name']] = self.parse_logs(f['f'])

    def parse_logs(self, f):
        data = {}
        for l in f.splitlines():
            s = l.split()
            data[s[0]] = s[1]
        return data
```

### Filtering by parsed sample names
MultiQC users can use the `--ignore-samples` flag to skip sample names
that match specific patterns. As sample names are generated in a different
way by every module, this filter has to be applied after log parsing.

There is a core function to do this task - assuming that your data is
in a dictionary with the first key as sample name, pass it through the
`self.ignore_samples` function as follows:

```python
self.yourdata = self.ignore_samples(self.yourdata)
```

This will remove any dictionary keys where the sample name matches
a user pattern.

### No files found
If your module cannot find any matching files, it needs to raise an
exception of type `UserWarning`. This tells the core MultiQC program
that no modules were found. For example:
```python
if len(self.mod_data) == 0:
    log.debug("Could not find any data in {}".format(config.analysis_dir))
    raise UserWarning
```

Note that this has to be raised as early as possible, so that it halts
the module progress. For example, if no logs are found then the module
should not create any files or try to do any computation.

### Custom sample names
Typically, sample names are taken from cleaned log filenames (the default
`f['s_name']` value returned). However, if possible, it's better to use
the name of the input file (allowing for concatenated log files).
To do this, you should use the `self.clean_s_name()` function, as
this will prepend the directory name if requested on the command line:

```python
input_fname = s[3] # Or parsed however
s_name = self.clean_s_name(input_fname, f['root'])
```

This function has already been applied to the contents of `f['s_name']`.

> `self.clean_s_name()` **must** be used on sample names parsed from the file
> contents. Without it, features such as prepending directories (`--dirs`)
> will not work.

### Identical sample names
If modules find samples with identical names, then the previous sample
is overwritten. It's good to print a log statement when this happens,
for debugging. However, most of the time it makes sense - programs often
create log files _and_ print to `stdout` for example.

```python
if f['s_name'] in self.bowtie_data:
    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
```

### Printing to the sources file
Finally, once you've found your file we want to add this information to the
`multiqc_sources.txt` file in the MultiQC report data directory. This lists
every sample name and the file from which this data came from. This is especially
useful if sample names are being overwritten as it lists the source used. This code
is typically written immediately after the above warning.

If you've used the `self.find_log_files` function, writing to the sources file
is as simple as passing the log file variable to the `self.add_data_source` function:
```python
for f in self.find_log_files('mymod'):
    self.add_data_source(f)
```

If you have different files for different sections of the module, or are
customising the sample name, you can tweak the fields. The default arguments
are as shown:
```python
self.add_data_source(f=None, s_name=None, source=None, module=None, section=None)
```

## Step 4 - Adding to the general statistics table
Now that you have your parsed data, you can start inserting it into the
MultiQC report. At the top of ever report is the 'General Statistics'
table. This contains metrics from all modules, allowing cross-module
comparison.

There is a helper function to add your data to this table. It can take
a lot of configuration options, but most have sensible defaults. At
it's simplest, it works as follows:

```python
data = {
    'sample_1': {
        'first_col': 91.4,
        'second_col': '78.2%'
    },
    'sample_2': {
        'first_col': 138.3,
        'second_col': '66.3%'
    }
}
self.general_stats_addcols(data, headers)
```

To give more informative table headers and configure things like
data scales and colour schemes, you can supply an extra dict:
```python
headers = OrderedDict()
headers['first_col'] = {
    'title': 'First',
    'description': 'My First Column',
    'scale': 'RdYlGn-rev'
}
headers['second_col'] = {
    'title': 'Second',
    'description': 'My Second Column',
    'max': 100,
    'min': 0,
    'scale': 'Blues',
    'suffix': '%'
}
self.general_stats_addcols(data, headers)
```

Here are all options for headers, with defaults:
```python
headers['name'] = {
    'namespace': '',                # Module name. Auto-generated for General Statistics.
    'title': '[ dict key ]',        # Short title, table column title
    'description': '[ dict key ]',  # Longer description, goes in mouse hover text
    'max': None,                    # Minimum value in range, for bar / colour coding
    'min': None,                    # Maximum value in range, for bar / colour coding
    'scale': 'GnBu',                # Colour scale for colour coding. Set to False to disable.
    'format': '{:,.1f}',            # Output format() string
    'shared_key': None              # See below for description
    'modify': None,                 # Lambda function to modify values
    'hidden': False                 # Set to True to hide the column on page load
}
```
* `namespace`
  * This prepends the column title in the mouse hover: _Namespace: Title_.
    It's automatically generated for the General Statistics table.
* `scale`
  * Colour scales are the names of ColorBrewer palettes. See below for available scales.
  * Add `-rev` to the name of a colour scale to reverse it
  * Set to `False` to disable colouring and background bars
* `shared_key`
  * Any string can be specified here, if other columns are found that share
    the same key, a consistent colour scheme and data scale will be used in
    the table. Typically this is set to things like `read_count`, so that
    the read count in a sample can be seen varying across analysis modules.
* `modify`
  * A python `lambda` function to change the data in some way when it is
    inserted into the table.
* `hidden`
  * Setting this to `True` will hide the column when the report loads. It can
    then be shown through the _Configure Columns_ modal in the report. This can
    be useful when data could be sometimes useful. For example, some modules
    show "percentage aligned" on page load but hide "number of reads aligned".

The typical use for the `modify` string is to divide large numbers such as read counts,
to make them easier to interpret. If handling read counts, there are three config variables
that should be used to allow users to change the multiplier: `read_count_multiplier`,
`read_count_prefix` and `read_count_desc`. For example:

```python
'title': '{} Reads'.format(config.read_count_prefix),
'description': 'Number of reads ({})'.format(config.read_count_desc),
'modify': lambda x: x * config.read_count_multiplier,
```

The colour scales are from [ColorBrewer2](http://colorbrewer2.org/) and are named as follows:
![color brewer](images/cbrewer_scales.png)

A third parameter can be passed to this function, `namespace`. This is usually
not needed - MultiQC automatically takes the name of the module that is calling
the function and uses this. However, sometimes it can be useful to overwrite this.

## Step 5 - Writing data to a file
In addition to printing data to the General Stats, MultiQC modules typically
also write to text-files to allow people to easily use the data in downstream
applications. This also gives the opportunity to output additional data that
may not be appropriate for the General Statistics table.

Again, there is a base class function to help you with this - just supply it
with a dictionary and a filename:

```python
data = {
    'sample_1': {
        'first_col': 91.4,
        'second_col': '78.2%'
    },
    'sample_2': {
        'first_col': 138.3,
        'second_col': '66.3%'
    }
}
self.write_data_file(data, 'multiqc_mymod')
```

If your output has a lot of columns, you can supply the additional
argument `sort_cols = True` to have the columns alphabetically sorted.

This function will also pay attention to the default / command line
supplied data format and behave accordingly. So the written file could
be a tab-separated file (default), `JSON` or `YAML`.

Note that any keys with more than 2 levels of nesting will be ignored
when being written to tab-separated files.

## Step 6 - Create report sections
Great! It's time to start creating sections of the report with more information.
To do this, use the `self.add_section()` helper function:

```python
self.add_section (
    name = 'First Module Section',
    anchor = 'mymod-first',
    description = 'My amazing module output, from the first section',
    help = "If you're not sure how to interpret the data, we can help!",
    plot = bargraph.plot(data)
)
self.add_section (
    name = 'Second Module Section',
    anchor = 'mymod-second',
    plot = linegraph.plot(data2)
)
self.add_section (
    content = '<p>Some custom HTML.</p>'
)
```
These will automatically be labelled and linked in the navigation (unless
the module has only one section or `name` is not specified).

## Step 7 - Plot some data
Ok, you have some data, now the fun bit - visualising it! Each of the plot
types is described in the _Plotting Functions_ section of the docs.

## Appendices

### Profiling Performance
It's important that MultiQC runs quickly and efficiently, especially on big
projects with large numbers of samples. The recommended method to check this is
by using `cProfile` to profile the code execution. To do this, run MultiQC as follows:

```bash
python -m cProfile -o multiqc_profile.prof /path/to/MultiQC/scripts/multiqc -f .
```

You can create a `.bashrc` alias to make this easier to run:
```bash
alias profile_multiqc='python -m cProfile -o multiqc_profile.prof /path/to/MultiQC/scripts/multiqc '
profile_multiqc -f .
```

MultiQC should run as normal, but produce the additional binary file `multiqc_profile.prof`.
This can then be visualised with software such as [SnakeViz](https://jiffyclub.github.io/snakeviz/).

To install SnakeViz and visualise the results, do the following:
```bash
pip install snakeviz
snakeviz multiqc_profile.prof
```

A web page should open where you can explore the execution times of different nested functions.
It's a good idea to run MultiQC with a comparable number of results from other tools (eg. FastQC)
to have a reference to compare against for how long the code should take to run.


### Adding Custom CSS / Javascript
If you would like module-specific CSS and / or JavaScript added to the template,
just add to the `self.css` and `self.js` dictionaries that come with the
`BaseMultiqcModule` class. The key should be the filename that you want your file to
have in the generated report folder _(this is ignored in the default template, which
includes the content file directly in the HTML)_. The dictionary value should be
the path to the desired file. For example, see how it's done in the FastQC module:

```python
self.css = {
    'assets/css/multiqc_fastqc.css' :
        os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_fastqc.css')
}
self.js = {
    'assets/js/multiqc_fastqc.js' :
        os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_fastqc.js')
}
```
