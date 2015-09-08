# Contributing to MultiQC
Suggestions to the developer are welcome! Contributions to the code are
even more welcome ;)

## Contribution workflow:
1. Create an issue describing the bug / suggestion / improvement / ... [here](https://github.com/ewels/MultiQC/issues).
    1. Specify that you're planning to work on this yourself.
    This avoids someone else duplicating your work by accident.
2. Fork this repository to your GitHub account
3. Make the necessary changes / additions within your forked repository
4. Submit a Pull Request and wait for the code to be reviewed and merged.

## How to get started
The basic prerequisite for a new module is a new directory in `/multiqc/`
with a `__init__.py` file. The directory should have the name of the
module.

For a starting place, I recommend that you have a look at the `star`
module, which is fairly simple. For more advanced examples see `fastq_screen`
(with a custom plotting function) or if you're feeling really brave,
`fastqc` (lots of custom CSS and JS).

Note - you can add other required files to your module directory. Note that
any subdirectories will need to contain an empty `__init__.py` file to be
copied when the package is installed. Your module will also need to copy
required files to the report at run time. See the `init_modfiles()` function
in the FastQC module for an example of this.

## General Statistics Table
The general statistics table at the top of every MultiQC report is a special
case section, designed to bring together numbers for samples from across
different modules.

To add to the table, update the `report` python object with the following keys:

```python
report['general_stats']['headers']['YOUR_KEY'] = '<th>My header</th>'
report['general_stats']['rows']['SAMPLE_NAME']['YOUR_KEY'] = '<td>Value</td>'
```

To customise display and colouring of your data, you can include certain
special HTML attributes to the content. The header `<th>` tag can have the
following attributes to apply colouring. The colours are specified in
`data-chroma-scale` and should be a scale specified by
[chroma.js](https://github.com/gka/chroma.js/wiki/Predefined-Colors)
`-rev` may be appended to the name to reverse it, eg. `RdYlGn-rev`
```html
class="chroma-col"
data-chroma-scale="SCALE"
data-chroma-max="100"
data-chroma-min="0"
```

It's best to keep the table header text short so that the column is narrow.
You can add mouse hover text with more detail using the following HTML:
```html
<th>
    <span data-toggle="tooltip" title="My Module: A longer description of what this column is">
        My Header
    </span>
</th>
```

Finally, you can style the data table cell as well, for instance by aligning
the text to the right:
```html
<td class="text-right">96.7%</td>
```

## Python Helper Functions
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

### self.dict_to_csv (data, delim="\t")
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

### self.plot_xy_data (data, config)
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

### self.plot_bargraph (data, cats, config)
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

Finally, a config dict can be supplied. The defaults are as follows:
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

## Javascript Helper Functions
The javascript bundled in the default MultiQC template has a number of
helper functions to make your life easier.

Note that the above Python functions make use of these, so
you typically won't need to call these directly.

### plot_xy_line_graph (target, data, config)

Plots a line graph with multiple series of (x,y) data pairs.

Available config options with default vars:
```javascript
config = {
    title: undefined,           // Plot title
    xlab: undefined,            // X axis label
    ylab: undefined,            // Y axis label
    xmax: undefined,            // Max x limit
    xmin: undefined,            // Min x limit
    xDecimals: true,            // Set to false to only show integer labels
    ymax: undefined,            // Max y limit
    ymin: undefined,            // Min y limit
    yDecimals: true,            // Set to false to only show integer labels
    tt_label: '{point.x}',      // Use to customise tooltip label, eg. '{point.x} base pairs'
    click_func: function(){},   // Javascript function to be called when a point is clicked
    cursor: undefined           // CSS mouse cursor type. Defaults to pointer when 'click_func' specified
}
```

An example of the markup expected, with the function being called:
```html
<div id="my_awesome_line_graph" class="hc-plot"></div>
<script type="text/javascript">
    my_data = [
        {
            name: 'Sample 1',
            data: [[1, 1.5], [1.5, 3.1], [2, 6.4]]
        },
        {
            name: 'Sample 2',
            data: [[1, 1.7], [1.5, 4.3], [2, 8.4]]
        },
    ];
    var my_config = {
        "title": "Best Plot Ever",
        "ylab": "Pings",
        "xlab": "Pongs"
    };
    $(function () {
        plot_xy_line_graph("#my_awesome_line_graph", my_data, my_config);
    });
</script>
```


### plot_stacked_bar_graph (target, names, data, config)

Plots a bar graph with multiple series containing multiple categories.

All available config options with default vars:
```javascript
config = {
    title: undefined,           // Plot title
    xlab: undefined,            // X axis label
    ylab: undefined,            // Y axis label
    ymax: undefined,            // Max y limit
    ymin: undefined,            // Min y limit
    yDecimals: true,            // Set to false to only show integer labels
    stacking: 'normal',         // Set to null to have category bars side by side (None in python)
    use_legend: true,           // Show / hide the legend
    click_func: undefined,      // Javascript function to be called when a point is clicked
    cursor: undefined,          // CSS mouse cursor type. Defaults to pointer when 'click_func' specified
}
```

An example of the markup expected, with the function being called:
```html
<div id="my_awesome_bar_plot" class="hc-plot"></div>
<script type="text/javascript">
    my_sample_names = ['Sample 1', 'Sample 2']
    my_sample_data = [{"data": [4, 7], "name": "Passed Test"}, {"data": [2, 3], "name": "Failed Test"}]
    var my_config = {
        "title": "My Awesome Plot",
        "ylab": "# Observations",
        "ymin": 0,
        "stacking": "normal"
    };
    $(function () {
        plot_stacked_bar_graph("#my_awesome_bar_plot", my_sample_names, my_sample_data, my_config);
    });
</script>
```

### Switching Between Counts and Percentages
If you're using the plotting functions above, it's easy to add a button which
switches between percentages and counts. Just add the following HTML above
your plot:
```html
<div class="btn-group switch_group">
    <button class="btn btn-default btn-sm active" data-action="set_numbers" data-target="#my_plot">Counts</button>
    <button class="btn btn-default btn-sm" data-action="set_percent" data-target="#my_plot">Percentages</button>
</div>
```
_NB:_ This markup is generated automatically with the Python `self.plot_bargraph()` function.



### Switching Plot Datasets
Much like the counts / percentages buttons above, you can add a button which
switches the data displayed in a single plot. Make sure that both datasets
are stored in named javascript variables, then add the following markup:
```html
<div class="btn-group switch_group">
    <button class="btn btn-default btn-sm active" data-action="set_data" data-ylab="First Data" data-newdata="data_var_1" data-target="#my_plot">Data 1</button>
    <button class="btn btn-default btn-sm" data-action="set_data" data-ylab="Second Data" data-newdata="data_var_2" data-target="#my_plot">Data 2</button>
</div>
```
Note the CSS class `active` which specifies which button is 'pressed' on page load.
`data-ylab` and `data-xlab` can be used to specify the new axes labels.
`data-newdata` should be the name of the javascript object with the new data
to be plotted and `data-target` should be the CSS selector of the plot to change.

### Custom Event Triggers
Some of the events that take place in the general javascript
code trigger jQuery events which you can hook into from within your
module's code. This allows you to take advantage of events generated
by the global theme whilst keeping your code modular.

```javascript
$(document).on('mqc_highlights', function(e, f_texts, f_cols, regex_mode){
    // This trigger is called when the highlight strings are
    // updated. Three variables are given - an array of search
    // strings (f_texts), an array of colours with corresponding
    // indexes (f_cols) and a boolean var saying whether the
    // search should be treated as a string or a regex (regex_mode)
});

$(document).on('mqc_hidesamples', function(e, f_texts, regex_mode){
    // This trigger is called when the Hide Samples filters change.
    // Two variables are given - an array of search strings
    // (f_texts) and a boolean saying whether the search should
    // be treated as a string or a regex (regex_mode)
});

$('#YOUR_PLOT_ID').on('mqc_plotresize', function(){
    // This trigger is called when a plot handle is pulled,
    // resizing the height
});
```

