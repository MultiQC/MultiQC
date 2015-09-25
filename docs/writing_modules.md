---
title: Writing Modules
layout: toc
---

# Writing Modules

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