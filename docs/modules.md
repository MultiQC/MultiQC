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
see the [plugins page](plugins.md)

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
        'qualimap = multiqc.modules.qualimap:MultiqcModule',
        [...]
    ]
}
```
Copy one of the existing module lines and change it to use your module name.
Once this is done, you'll need to update your installation of MultiQC:
```
python setup.py develop
```

### MultiQC config
The final place that needs to know about your module is the MultiQC config.
This is required as MultiQC modules need to be run in the correct order.

In `multiqc/utils/config.py` you should see a list variable called `module_order`.
This contains the name of modules in order of precedence. Add your module here
in an appropriate position.

### MultiqcModule Class
If you've copied one of the other entry point statements, it will have
ended in `:MultiqcModule` - this tells MultiQC to try to execute a class or
function called `MultiqcModule`.

To use the helper functions bundled with MultiQC, you should extend this
class from `multiqc.BaseMultiqcModule`. This will give you access to a
number of functions on the `self` namespace. For example:
```python
class MultiqcModule(multiqc.BaseMultiqcModule):
    def __init__(self, report):
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

You should see the patterns for all other modules to give you an idea,
but you want a yaml key with the name of your module, then either `fn`
or `contents` for strings to match against filenames or file contents:
```yaml
mymod:
    fn: _myprogram.txt
myothermod:
    contents: This is myprogram v1.3
```

Note that if you want to find multiple log files, you can nest these
dictionaries (though they must end with either `fn` or `contents`).
For example, see the `FastQC` module:
```yaml
fastqc:
    data:
        fn: fastqc_data.txt
    zip:
        fn: _fastqc.zip
```

You can supply a list of strings if needed, eg. the `bismark` module:
```yaml
bismark:
    align:
        fn:
            - _PE_report.txt
            - _SE_report.txt
```

The value of adding these strings here is that they can be overwritten
by users in their own config files. This is helpful as people have weird
and wonderful processing pipelines with their own file naming conventions.

Once your strings are added, you can call them in your module with the
`config.sp['mymod']`. Next, use the base function `self.find_log_files()`
to look for your files like this:
```python
self.find_log_files(config.sp['mymod'], filehandles=False)
```

This will recursively search the analysis directories looking for a matching
file name (if the `fn` key is there) or a text string held within a file
(if the `contents` key is there). Contents matching is only done on files
smaller than `config.log_filesize_limit` (default 1MB).
Note that both `fn` and `contents` can be used in combination
if required - files will be returned if anything matches (`OR` not `AND`).

This function yields a dictionary with various information about matching
files. The `f` key contains the contents of matching files by default.
```python
# Find all files for mymod
for f in self.find_log_files(config.sp['mymod']):
    print f['f']        # File contents
    print f['s_name']   # Sample name (from cleaned fn)
    print f['root']     # Directory file was in
    print f['fn']       # Filename
```

If `filehandles=True` is specified, the `f` key contains a file handle
instead:
```python
# Find all files which contain the string 'My Statistic:'
# Return a filehandle instead of the file contents
for f in self.find_log_files(config.sp['mymod'], filehandles=True):
    line = f['f'].readline()  # f['f'] is now a filehandle instead of contents
```

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
        for f in self.find_log_files(config.sp['mymod']):
            self.mod_data[f['s_name']] = self.parse_logs(f['f'])
    
    def parse_logs(self, f):
        data = {}
        for l in f.splitlines():
            s = l.split()
            data[s[0]] = s[1]
        return data
```

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
for f in self.find_log_files(config.sp['mymod']):
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
    'format': '{:.1f}%'
}
self.general_stats_addcols(data, headers)
```

Here are all options for headers, with defaults:
```python
headers['name'] = {
    'title': '[ dict key ]',        # Short title, table column title
    'description': '[ dict key ]',  # Longer description, goes in mouse hover text
    'max': None,                    # Minimum value in range, for bar / colour coding
    'min': None,                    # Maximum value in range, for bar / colour coding
    'scale': 'GnBu',                # Colour scale for colour coding
    'format': '{:.1f}',             # Output format() string
    'shared_key': None              # See below for description
    'modify': None,                 # Lambda function to modify values
}
```
* `scale`
  * Colour scales are the names of ColorBrewer palettes. See the
    [chroma.js documentation](https://github.com/gka/chroma.js/wiki/Predefined-Colors)
    for a list of available colour scales
  * Add `-rev` to the name of a colour scale to reverse it
* `shared_key`
  * Any string can be specified here, if other columns are found that share
    the same key, a consistent colour scheme and data scale will be used in
    the table. Typically this is set to things like `read_count`, so that
    the read count in a sample can be seen varying across analysis modules.
* `modify`
  * A python `lambda` function to change the data in some way when it is
    inserted into the table. Typically, this is used to divide numbers to
    show millions: `'modify': lambda x: x / 1000000`
Example of supplied data and headers dicts:

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
If you only have one plot / section to create, just add it to the introduction.
For example (content to be replaced with a brilliant plot in the bit of the docs):

```python
self.intro += 'My amazing module output'
```

If you have multiple plots to show (_eg._ the Qualimap and FastQC modules),
you can create a list to hold the sections:

```python
self.sections = list()
self.sections.append({
    'name': 'First Module Section',
    'anchor': 'mymod_first',
    'content': 'My amazing module output, from the first section'
})
self.sections.append({
    'name': 'Second Module Section',
    'anchor': 'mymod_second',
    'content': 'My amazing module output, from the second section'
})
```
The will automatically be labelled and linked in the navigation.

## Step 7 - Plotting bar graphs
Simple data can be plotted in bar graphs. Many MultiQC modules make use
of stacked bar graphs. Here, the `self.plot_bargraph()` module comes to
the rescue. A basic example is as follows:
```python
data = {
    'sample 1': {
        'aligned': 23542,
        'not_aligned': 343,
    },
    'sample 2': {
        'not_aligned': 7328,
        'aligned': 1275,
    }
}
html_content = self.plot_bargraph(data)
```

To specify the order of categories in the plot, you can supply a list of
dictionary keys. This can also be used to exclude a key from the plot.

```python
cats = ['aligned', 'not_aligned']
html_content = self.plot_bargraph(data, cats)
```

If `cats` is given as a dict instead of a list, you can specify a nice name
and a colour too. Make it an OrderedDict to specify the order:
```python
from collections import OrderedDict
cats = OrderedDict()
cats['aligned'] = {
    'name': 'Aligned Reads',
    'color': '#8bbc21'
}
cats['not_aligned'] = {
    'name': 'Unaligned Reads',
    'color': '#f7a35c'
}
```

Finally, a third variable can be supplied with configuration variables for
the plot. The defaults are as follows:
```python
config = {
    # Building the plot
    'id': '<random string>',                # HTML ID used for plot
    'cpswitch': True,                       # Show the 'Counts / Percentages' switch?
    'cpswitch_c_active': True,              # Initial display with 'Counts' specified? False for percentages.
    'cpswitch_counts_label': 'Counts',      # Label for 'Counts' button
    'cpswitch_percent_label': 'Percentages' # Label for 'Percentages' button
    'logswitch': False,                     # Show the 'Log10' switch?
    'logswitch_active': False,              # Initial display with 'Log10' active?
    'logswitch_label': 'Log10',             # Label for 'Log10' button
    # Customising the plot
    'title': None,                          # Plot title
    'xlab': None,                           # X axis label
    'ylab': None,                           # Y axis label
    'ymax': None,                           # Max y limit
    'ymin': None,                           # Min y limit
    'yDecimals': True,                      # Set to false to only show integer labels
    'ylab_format': None,                    # Format string for x axis labels. Defaults to {value}
    'stacking': 'normal',                   # Set to None to have category bars side by side
    'use_legend': True,                     # Show / hide the legend
    'click_func': None,                     # Javascript function to be called when a point is clicked
    'cursor': None,                         # CSS mouse cursor type.
    'tt_percentages': True,                 # Show the percentages of each count in the tooltip
}
```

### Switching datasets
It's possible to have single plot with buttons to switch between different
datasets. To do this, give a list of data objects (same formats as described
above). Also add the following config options to supply names to the buttons
and graph labels:
```python
config = {
    'data_labels': [
        {'name': 'Reads', 'ylab': 'Number of Reads'},
        {'name': 'Frags', 'ylab': 'Percentage of Fragments', 'ymax':100}
    ]
}
```
If supplying multiple datasets, you can also supply a list of category
objects. Make sure that they are in the same order as the data. If not
supplied, these will be guessed from the data keys. See the bismark module
plots for an example of this in action.

### Interactive / Flat image plots
Note that the `self.plot_bargraph()` function can generate both interactive
JavaScript (HighCharts) powered report plots _and_ flat image plots made using
MatPlotLib. This choice is made within the function based on config variables
such as number of dataseries and command line flags.

Note that both plot types should come out looking pretty much identical. If
you spot something that's missing in the flat image plots, let me know.


## Step 8 - Plotting line graphs
This base function works much like the above, but for two-dimensional
data, to produce line graphs. It expects a dictionary in the following format:
```python
data = {
    'sample 1': {
        '<x val 1>': '<y val 1>',
        '<x val 2>': '<y val 2>',
    },
    'sample 2': {
        '<x val 1>': '<y val 1>',
        '<x val 2>': '<y val 2>',
    }
}
html_content = self.plot_xy_data(data)
```

Additionally, a config dict can be supplied. The defaults are as follows:
```python
config = {
    # Building the plot
    'id': '<random string>',     # HTML ID used for plot    
    'categories': <anything>,    # Set key to use x values as categories instead of numbers.
    'colors': dict()             # Provide dict with keys = sample names and values colours
    'extra_series': None,        # See section below
    # Plot configuration
    'title': None,               # Plot title
    'xlab': None,                # X axis label
    'ylab': None,                # Y axis label
    'xCeiling': None,            # Maximum value for automatic axis limit (good for percentages)
    'xFloor': None,              # Minimum value for automatic axis limit
    'xMinRange': None,           # Minimum range for axis
    'xmax': None,                # Max x limit
    'xmin': None,                # Min x limit
    'xDecimals': True,           # Set to false to only show integer labels
    'yCeiling': None,            # Maximum value for automatic axis limit (good for percentages)
    'yFloor': None,              # Minimum value for automatic axis limit
    'yMinRange': None,           # Minimum range for axis
    'ymax': None,                # Max y limit
    'ymin': None,                # Min y limit
    'yDecimals': True,           # Set to false to only show integer labels
    'yPlotBands': None,          # Highlighted background bands. See http://api.highcharts.com/highcharts#yAxis.plotBands
    'xPlotBands': None,          # Highlighted background bands. See http://api.highcharts.com/highcharts#xAxis.plotBands
    'yPlotLines': None,          # Highlighted background bands. See http://api.highcharts.com/highcharts#yAxis.plotLines
    'xPlotLines': None,          # Highlighted background bands. See http://api.highcharts.com/highcharts#xAxis.plotLines
    'tt_label': '{point.x}: {point.y:.2f}', # Use to customise tooltip label, eg. '{point.x} base pairs'
    'pointFormat': None,         # Replace the default HTML for the entire tooltip label
    'click_func': function(){},  # Javascript function to be called when a point is clicked
    'cursor': None               # CSS mouse cursor type. Defaults to pointer when 'click_func' specified
    'reversedStacks': False      # Reverse the order of the category stacks. Defaults True for plots with Log10 option
}
html_content = self.plot_xy_data(data, config)
```

### Switching datasets
You can also have a single plot with buttons to switch between different
datasets. To do this, just supply a list of data dicts instead (same
formats as described above). Also add the following config options to
supply names to the buttons and graph labels:
```python
config = {
    'data_labels': [
        {'name': 'DS 1', 'ylab': 'Dataset 1'},
        {'name': 'DS 2', 'ylab': 'Dataset 2'}
    ]
}
```
All of these config values are optional, the function will default
to sensible values if things are missing. See the cutadapt module
plots for an example of this in action.

### Additional data series
Sometimes, it's good to be able to specify specific data series manually.
To do so, set `config['extra_series']` as a `list` of `dict`s. For example,
the Preseq module does this to create the dotted `x = y` reference line:
```python
config = {
    'extra_series': [
        {
            'name': 'x = y',
            'data': [[0, 0], [maxval, maxval]],
            'dashStyle': 'Dash',
            'lineWidth': 1,
            'color': '#000000',
            'marker': { 'enabled': False },
            'enableMouseTracking': False,
            'showInLegend': False,
        }
    ]
}
html_content = self.plot_xy_data(data, config)
``` 

## Appendix A - Including module-specific files
Occasionally, a module may require additional files to be copied at
the time of report generation. To do this, use the `copy_module_files`
function. `files` is a list of file paths, the second argument is the
path that the files should be relative to (usually `__file__`). For example:
```python
self.copy_module_files(['includes/file1.txt', 'includes/file2.txt'], __file__)
```

### Appendix B - Adding Custom CSS / Javascript
If you would like module-specific CSS and / or JavaScript added to the template,
create lists in the module called `self.css` and `self.js` with the relative file
paths - these will be looped over in the template file and the files included
appropriately. Make sure that the files themselves are copied using the above
method. The FastQC module does this to generate it's custom plots:

```python
self.css = [ os.path.join('assets', 'css', 'multiqc_fastqc.css') ]
self.js = [ os.path.join('assets', 'js', 'multiqc_fastqc.js') ]
self.copy_module_files(self.css + self.js, __file__)
```

### Appendix C - Custom plotting functions
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

These plotting functions aren't very helpful, but hopefully you get the idea.
