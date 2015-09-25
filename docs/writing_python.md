---
title: Python Help
layout: toc
---

# Python Helper Functions
When defining your new Python module, you should inherit the base module
class as follows:

```python
class MultiqcModule(multiqc.BaseMultiqcModule):
    def __init__(self, report):
        # Initialise the parent object
        super(MultiqcModule, self).__init__()
```

Doing so gives you access to common Python functions. These are currently
few in number (just one) but may be extended in the near future.

## Writing data to a file

> `self.dict_to_csv (data, delim="\t")`

This function takes a 2D dictionary and returns a string suitable for
writing to a .csv file. First key should be sample name (row header),
second key should be field (column header).

The function takes a dictionary as input, plus an optional 'delim'
field to specify column delimiter (default: tab).

You can see an example of this function in the featureCounts module:
```python
# Write parsed report data to a file
with open (os.path.join(self.output_dir, 'report_data', 'multiqc_featureCounts.txt'), "w") as f:
    print( self.dict_to_csv( self.featurecounts_data ), file=f)
```

## Plotting line graphs

> `self.plot_xy_data (data, config, original_plots)`

This function takes a dict of data and plots an XY line graph.
It expects a dictionary in the following format:
```python
data = {
    'sample 1': {
        '<x val 1>': '<y val 1>',
        '<x val 2>': '<y val 2>',
        [ .. ]
    },
    'sample 2': {
        '<x val 1>': '<y val 1>',
        '<x val 2>': '<y val 2>',
        [ .. ]
    },
    [ .. ]
}
html += self.plot_xy_data(data)
```
Additionally, a config dict can be supplied. The defaults are as follows:
```python
config = {
    'id': '<random string>', # HTML ID used for plot    
}
```
This dictionary can also have any of the javascript config options.
See the [plot_xy_line_graph()](CONTRIBUTING.md#plot_xy_line_graph-target-data-config) below for those.

### Switching datasets
You can also have a single plot with buttons to switch between different
datasets. To do this, just supply a list of data dicts instead (same
formats as described above). Also add the following config options to
supply names to the buttons and graph labels:
```python
config = {
    'data_labels': [{'name': 'DS 1', 'ylab': 'Dataset 1'},
                    {'name': 'DS 2', 'ylab': 'Dataset 2'}] 
}
```
All of these config values are optional, the function will default
to sensible values if things are missing. See the cutadapt module
plots for an example of this in action.

### Showing original plots
Finally, it's possible to get the plots to show original data sources
when clicked (for an example, see the FastQC plots). To achieve this
magic, make sure that the images are copied to the report and pass
the `original_plots` variable to the function. It should look like this:
```python
original_plots = [
    {
        's_name': "My First Sample",
        'img_path': 'report_data/my_module/first_sample_original_image.png'
    },
    {
        's_name': "My Second Sample",
        'img_path': 'report_data/my_module/second_sample_original_image.png'
    },
]
```
Note that the _Prev_ / _Next_ buttons will cycle through these images in the
order supplied, so it makes sense to sort by sample name alphabetically.

## Plotting bar graphs

> `self.plot_bargraph (data, cats, config)`

Takes a dict of data and plots a bar graph. The expected data structure
is as follows:
```python
data = {
    'sample 1': {
        'category 1': '<val>',
        'category 2': '<val>',
        [ .. ]
    },
    'sample 2': {
        'category 1': '<val>',
        'category 2': '<val>',
        [ .. ]
    },
    [ .. ]
}
html += self.plot_bargraph(data)
```
Cats is an optional extra which allows you to specify how the different
categories are displayed. If not specified, the categories are determined
from the sample keys, but their order is random.

If given as a list you can dictates the order, and ignore certain keys:
```python
cats = ['first_category', 'second_category', 'fourth_category']
```

If given as a dict, you can specify a nice name and a colour (optional).
Make it an OrderedDict to specify the order as well:
```python
from collections import OrderedDict
cats = OrderedDict()
cats['first_category'] = {
    'name': 'The First Category',
    'color': '#8bbc21'
}
cats['second_category'] = {
    'name': 'The Second Category',
    'color': '#f7a35c'
}
```

Like the X,Y plot, a config dict can be supplied. The defaults are as follows:
```python
config = {
    'id': '<random string>',                # HTML ID used for plot
    'cpswitch': True,                       # Show the 'Counts / Percentages' switch?
    'cpswitch_c_active': True,              # Initial display with 'Counts' specified?
    'cpswitch_counts_label': 'Counts',      # Label for 'Counts' button
    'cpswitch_percent_label': 'Percentages' # Label for 'Percentages' button
}
```
This dictionary can also have any of the javascript config options.
See the [plot_stacked_bar_graph()](CONTRIBUTING.md#plot_stacked_bar_graph-target-names-data-config)
docs below for those.

### Switching datasets
You can also have a single plot with buttons to switch between different
datasets. To do this, give a list of data objects (same formats as described
above). Also add the following config options to supply names to the buttons
and graph labels:
```python
config = {
    'data_labels': [{'name': 'DS 1', 'ylab': 'Dataset 1'},
                    {'name': 'DS 2', 'ylab': 'Dataset 2'}] 
}
```
If supplying multiple datasets, you can also supply multiple category
objects. Make sure that they are in the same order as the data. If not
supplied, these will be guessed from the data keys.

See the bismark module plots for an example of this in action.