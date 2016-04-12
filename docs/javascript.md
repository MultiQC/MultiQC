# Javascript Functions
The javascript bundled in the default MultiQC template has a number of
helper functions to make your life easier.

Note that the base module Python functions make use of these, so
you typically won't need to call these directly.

## Plotting line graphs
`plot_xy_line_graph (target, ds)`

Plots a line graph with multiple series of (x,y) data pairs. Used by
the [self.plot_xy_data()](CONTRIBUTING.md#selfplot_xy_data-data-config)
python function.

Data and configuration must be added to the document level
`mqc_plots` variable on page load, using the target as the key.
The variables used are as follows:
```javascript
mqc_plots[target]['plot_type'] = 'xy_line';
mqc_plots[target]['config'];
mqc_plots[target]['datasets'];
```

Multiple datasets can be added in the `['datasets']` array. The supplied
variable `ds` specifies which is plotted (defaults to `0`).

Available config options with default vars:
```javascript
config = {
    title: undefined,            // Plot title
    xlab: undefined,             // X axis label
    ylab: undefined,             // Y axis label
    xCeiling: undefined,         // Maximum value for automatic axis limit (good for percentages)
    xFloor: undefined,           // Minimum value for automatic axis limit
    xMinRange: undefined,        // Minimum range for axis
    xmax: undefined,             // Max x limit
    xmin: undefined,             // Min x limit
    xDecimals: true,             // Set to false to only show integer labels
    yCeiling: undefined,         // Maximum value for automatic axis limit (good for percentages)
    yFloor: undefined,           // Minimum value for automatic axis limit
    yMinRange: undefined,        // Minimum range for axis
    ymax: undefined,             // Max y limit
    ymin: undefined,             // Min y limit
    yDecimals: true,             // Set to false to only show integer labels
    yPlotBands: undefined,       // Highlighted background bands. See http://api.highcharts.com/highcharts#yAxis.plotBands
    xPlotBands: undefined,       // Highlighted background bands. See http://api.highcharts.com/highcharts#xAxis.plotBands
    tt_label: '{point.x}: {point.y:.2f}', // Use to customise tooltip label, eg. '{point.x} base pairs'
    pointFormat: undefined,      // Replace the default HTML for the entire tooltip label
    click_func: function(){},    // Javascript function to be called when a point is clicked
    cursor: undefined            // CSS mouse cursor type. Defaults to pointer when 'click_func' specified
}
```

An example of the markup expected, with the function being called:
```html
<div id="my_awesome_line_graph" class="hc-plot"></div>
<script type="text/javascript">
    mqc_plots['#my_awesome_bar_plot']['plot_type'] = 'xy_line';
    mqc_plots['#my_awesome_line_graph']['datasets'] = [
        {
            name: 'Sample 1',
            data: [[1, 1.5], [1.5, 3.1], [2, 6.4]]
        },
        {
            name: 'Sample 2',
            data: [[1, 1.7], [1.5, 4.3], [2, 8.4]]
        },
    ];
    mqc_plots['#my_awesome_line_graph']['config'] = {
        "title": "Best Plot Ever",
        "ylab": "Pings",
        "xlab": "Pongs"
    };
    $(function () {
        plot_xy_line_graph('#my_awesome_line_graph');
    });
</script>
```

## Plotting bar graphs
`plot_stacked_bar_graph (target, ds)`

Plots a bar graph with multiple series containing multiple categories.
Used by the [self.plot_bargraph()](CONTRIBUTING.md#selfplot_bargraph-data-cats-config)
python function.

Data and configuration must be added to the document level
`mqc_plots` variable on page load, using the target as the key.
The variables used are as follows:
```javascript
mqc_plots[target]['plot_type'] = 'bar_graph';
mqc_plots[target]['config'];
mqc_plots[target]['datasets'];
mqc_plots[target]['samples'];
```

All available config options with default vars:
```javascript
config = {
    title: undefined,           // Plot title
    xlab: undefined,            // X axis label
    ylab: undefined,            // Y axis label
    ymax: undefined,            // Max y limit
    ymin: undefined,            // Min y limit
    yDecimals: true,            // Set to false to only show integer labels
    ylab_format: undefined,     // Format string for x axis labels. Defaults to {value}
    stacking: 'normal',         // Set to null to have category bars side by side (None in python)
    xtype: 'linear',            // Axis type. 'linear' or 'logarithmic'
    use_legend: true,           // Show / hide the legend
    click_func: undefined,      // Javascript function to be called when a point is clicked
    cursor: undefined,          // CSS mouse cursor type. Defaults to pointer when 'click_func' specified
    tt_percentages: true,       // Show the percentages of each count in the tooltip
    reversedStacks: false,      // Reverse the order of the categories in the stack.
}
```

An example of the markup expected, with the function being called:
```html
<div id="my_awesome_bar_plot" class="hc-plot"></div>
<script type="text/javascript">
    mqc_plots['#my_awesome_bar_plot']['plot_type'] = 'bar_graph';
    mqc_plots['#my_awesome_bar_plot']['samples'] = ['Sample 1', 'Sample 2']
    mqc_plots['#my_awesome_bar_plot']['datasets'] = [{"data": [4, 7], "name": "Passed Test"}, {"data": [2, 3], "name": "Failed Test"}]
    mqc_plots['#my_awesome_bar_plot']['config'] = {
        "title": "My Awesome Plot",
        "ylab": "# Observations",
        "ymin": 0,
        "stacking": "normal"
    };
    $(function () {
        plot_stacked_bar_graph("#my_awesome_bar_plot");
    });
</script>
```

## Switching counts and percentages
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



## Switching plot datasets
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

## Custom event triggers
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

$(document).on('mqc_renamesamples', function(e, f_texts, t_texts, regex_mode){
    // This trigger is called when samples are renamed
    // Three variables are given - an array of search
    // strings (f_texts), an array of replacements with corresponding
    // indexes (t_texts) and a boolean var saying whether the
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

$('#YOUR_PLOT_ID').on('mqc_original_series_click', function(e, name){
    // A plot able to show original images has had a point clicked.
    // 'name' contains the name of the series that was clicked
});

$('#YOUR_PLOT_ID').on('mqc_original_chg_source', function(e, name){
    // A plot with original images has had a request to change the
    // original image source (eg. pressing Prev / Next)
});
```