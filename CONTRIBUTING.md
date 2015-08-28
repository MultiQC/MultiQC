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

# Javascript Helper Functions
The javascript bundled in the default MultiQC template has a number of
helper functions to make your life easier.

The following minimal examples make use of these helper functions to
create plots which can be resized / have highlighted sample names and so on:

### Example Line Graph
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

All available config options with default vars:
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
    use_legend: true,           // Show / hide the legend
    click_func: function(){}    // Javascript function to be called when a point is clicked
}
```


### Example Bar Plot HTML
```html
<div id="my_awesome_bar_plot" class="hc-plot"></div>
<script type="text/javascript"> \n\
    my_sample_names = ['sample 1', 'sample 2']
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

All available config options with default vars:
```javascript
config = {
    title: undefined,           // Plot title
    xlab: undefined,            // X axis label
    ylab: undefined,            // Y axis label
    ymax: undefined,            // Max y limit
    ymin: undefined,            // Min y limit
    yDecimals: true,            // Set to false to only show integer labels
    stacking: undefined,        // Set to 'normal' or 'percent' to stack http://api.highcharts.com/highcharts#series<bar>.stacking
    use_legend: true,           // Show / hide the legend
    tt_label: undefined         // Use to customise tooltip label, eg. '{point.x} base pairs'
}
```

### Custom Event Triggers
Some of the events that take place in the general javascript javascript
code trigger jQuery events which you can hook into from within your
module's code. This allows you to take advantage of events generated
by the global theme whilst keeping your code modular.

```javascript
$(document).on('mqc_highlights:reset', function(){
    // This trigger is called when something changes with the
    // sample name highlights. This could include a highlight
    // being removed, so should be used to reset plots
});
$(document).on('mqc_highlights:apply', function(f_text, f_col){
    // This trigger is called once for each highlight string
    // when something changes. Two variables are given - the
    // search string (f_text) and the colour (f_col)
});

$('#YOUR_PLOT_ID').on('mqc_plotresize', function(){
    // This trigger is called when a plot handle is pulled,
    // resizing the height
});
```



