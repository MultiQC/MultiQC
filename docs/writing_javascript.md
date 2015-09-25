---
title: Javascript Help
layout: toc
---

# Javascript Helper Functions
The javascript bundled in the default MultiQC template has a number of
helper functions to make your life easier.

Note that the above Python functions make use of these, so
you typically won't need to call these directly.

## Plotting line graphs

> `plot_xy_line_graph (target, data, config)`

Plots a line graph with multiple series of (x,y) data pairs. Used by
the [self.plot_xy_data()](CONTRIBUTING.md#selfplot_xy_data-data-config)
python function.

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
    tt_label: '{point.x}: {point.y:.2f}', // Use to customise tooltip label, eg. '{point.x} base pairs'
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

## Plotting bar graphs

> `plot_stacked_bar_graph (target, names, data, config)`

Plots a bar graph with multiple series containing multiple categories.
Used by the [self.plot_bargraph()](CONTRIBUTING.md#selfplot_bargraph-data-cats-config)
python function.

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

## Switching Between Counts and Percentages
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



## Switching Plot Datasets
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

## Custom Event Triggers
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

$('#YOUR_PLOT_ID').on('mqc_original_series_click', function(e, name){
    // A plot able to show original images has had a point clicked.
    // 'name' contains the name of the series that was clicked
});

$('#YOUR_PLOT_ID').on('mqc_original_chg_source', function(e, name){
    // A plot with original images has had a request to change the
    // original image source (eg. pressing Prev / Next)
});
```