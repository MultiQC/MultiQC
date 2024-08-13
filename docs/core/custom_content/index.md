---
title: Custom Content
description: Report on your data, even without a MultiQC module
---

# Introduction

Bioinformatics projects often include non-standardised analyses, with results from custom
scripts or in-house packages. It can be frustrating to have a MultiQC report describing
results from 90% of your pipeline but missing the final key plot. To help with this,
MultiQC has a special _"custom content"_ module.

Custom content parsing is a little more restricted than standard modules. Specifically:

- Only one plot per section is possible
- Plot customisation is more limited

All plot types can be generated using custom content - see the
[test files](https://github.com/MultiQC/test-data/tree/main/data/custom_content)
for examples of how data should be structured.

:::note
Use the name `custom_content` to refer to this module within configuration
settings that require a module name, such as [`module_order`](../reports/customisation.md#order-of-modules) or
[`run_modules`](../reports/customisation.md#removing-modules-or-sections).
:::

## Data from a released tool

If your data comes from a released bioinformatics tool, you shouldn't be using this
feature of MultiQC! Sure, you can probably get it to work, but it's better if a
fully-fledged core MultiQC module is written instead. That way, other users of MultiQC
can also benefit from results parsing.

Note that proper MultiQC modules are more robust and powerful than this custom-content
feature. You can also [write modules](../development/modules.md)
in [MultiQC plugins](../development/plugins.md) if they're not suitable for
general release.

## Images

As of MultiQC v1.7, you can import custom images into your MultiQC reports.
Simply add `_mqc` to the end of the filename for `.png`, `.jpg` or `.jpeg` files, for example:
`my_image_file_mqc.png` or `summmary_diagram.jpeg`.

Images will be embedded within the HTML file, so will be self contained.
Note that this means that it's very possible to make the HTML file very very large if abused!

The report section name and description will be automatically based on the filename.

Note that if you are using `sp:` to take in images with a custom filename you need to also set `ignore_images: false` in your config. For example:

```yaml
custom_data:
  my_custom_content_image:
    section_name: "My nice image"
sp:
  my_custom_content_image:
    fn: "*.png"
ignore_images: false
```

## MultiQC-specific data file

If you can choose exactly how your data output looks, then the easiest way to parse it
is to use a MultiQC-specific format. If the filename ends in `*_mqc.(yaml|yml|json|txt|csv|tsv|log|out|png|jpg|jpeg|html)`
then it will be found by any standard MultiQC installation with no additional customisation
required (v0.9 onwards).

These files contain configuration information specifying how the data should be parsed,
alongside the data. If you want to use YAML, this is an example of how it should look:

```yaml
id: "my_pca_section"
section_name: "PCA Analysis"
description: "This plot shows the first two components from a principal component analysis."
plot_type: "scatter"
pconfig:
  id: "pca_scatter_plot"
  title: "PCA Plot"
  xlab: "PC1"
  ylab: "PC2"
data:
  sample_1: { x: 12, y: 14 }
  sample_2: { x: 8, y: 6 }
  sample_3: { x: 5, y: 11 }
  sample_4: { x: 9, y: 12 }
```

:::note
This example YAML file is data only, and is not to be confused with a config file (though the two look very similar).
See the docs [Data as part of MultiQC config](#data-as-part-of-multiqc-config) for more on that.
:::

The file format can also be JSON:

```json
{
  "id": "custom_data_lineplot",
  "section_name": "Custom JSON File",
  "description": "This plot is a self-contained JSON file.",
  "plot_type": "linegraph",
  "pconfig": {
    "id": "custom_data_linegraph",
    "title": "Output from my JSON file",
    "ylab": "Number of things",
    "xDecimals": false
  },
  "data": {
    "sample_1": { "1": 12, "2": 14, "3": 10, "4": 7, "5": 16 },
    "sample_2": { "1": 9, "2": 11, "3": 15, "4": 18, "5": 21 }
  }
}
```

Note that if you're using `plot_type: html` then `data` just takes a string, with no sample keys.

For maximum compatibility with other tools, you can also use comma-separated or tab-separated files.
Include commented header lines with plot configuration in YAML format:

```bash
# id: "Output from my script'
# section_name: 'Custom data file'
# description: 'This output is described in the file header. Any MultiQC installation will understand it without prior configuration.'
# format: 'tsv'
# plot_type: 'bargraph'
# pconfig:
#    id: 'custom_bargraph_w_header'
#    ylab: 'Number of things'
Category_1    374
Category_2    229
Category_3    39
Category_4    253
```

You can easily inject custom HTML snippets by ending the filename with `_mqc.html` - again the
embedded config works in a similar way, but with a HTML comment:

```html
<!--
id: 'custom-html'
section_name: 'Custom HTML'
description: 'This section is created using a custom HTML file'
-->
<p>Some custom HTML content here.</p>
```

If no configuration is given, MultiQC will do its best to guess how to visualise your data appropriately.
To see examples of typical file structures which are understood, see the
[test data](https://github.com/MultiQC/test-data/tree/main/data/custom_content/no_config)
used to develop this code. Something will be probably be shown, but it may produce unexpected results.

:::note
Check [Tricky extras](#tricky-extras) for certain caveats about formatting headers for custom
`tsv` or `csv` files, particularly for the first column.
:::

## Data as part of MultiQC config

If you are already using a MultiQC config file to add data to your report (for example,
[titles / introductory text](../getting_started/config.md)), you can
give data within this file too. This can be in any MultiQC config file (for example,
passed on the command line with `-c my_yaml_file.yaml` or in your launch directory as
`multiqc_config.yml` - see [Configuration](../getting_started/config.md)).

:::note
This is not to be confused with the YAML data files described in the above section,
[MultiQC-specific data file](#multiqc-specific-data-file).
For example, MultiQC config files will _not_ be found with `_mqc.yml` file extensions.
:::

This is useful as you can
keep everything contained within a single file (including stuff unrelated to this
specific _custom content_ feature of MultiQC).

To be understood by MultiQC, the `custom_data` key must be found.
This must contain a section with a unique id, specific to your new report section.
Finally, the contents of this second dictionary will look the same as the above
stand-alone `YAML` files. For example:

```yaml
custom_data:
  my_data_type:
    id: "mqc_config_file_section"
    section_name: "My Custom Section"
    description: "This data comes from a single multiqc_config.yaml file"
    plot_type: "bargraph"
    pconfig:
      id: "barplot_config_only"
      title: "MultiQC Config Data Plot"
      ylab: "Number of things"
    data:
      sample_a:
        first_thing: 12
        second_thing: 14
      sample_b:
        first_thing: 8
        second_thing: 6
      sample_c:
        first_thing: 11
        second_thing: 5
      sample_d:
        first_thing: 12
        second_thing: 9
```

Or to add data to the General Statistics table:

```yaml
custom_data:
  my_genstats:
    plot_type: "generalstats"
    pconfig:
      - col_1:
          max: 100
          min: 0
          scale: "RdYlGn"
          suffix: "%"
      - col_2:
          min: 0
    data:
      sample_a:
        col_1: 14.32
        col_2: 1.2
      sample_b:
        col_1: 84.84
        col_2: 1.9
```

:::note
Use a **list** of headers in `pconfig` (keys prepended with `-`) to specify the order
of columns in the General Statistics table.
:::

See the [general statistics docs](../development/modules.md#step-3---adding-to-the-general-statistics-table)
for more information about configuring data for the General Statistics table.

## Separate configuration and data files

It's not always possible or desirable to include MultiQC configuration within a data file.
If this is the case, you can add to the MultiQC configuration to specify how input files
should be parsed.

As described in the [Data as part of MultiQC config](#data-as-part-of-multiqc-config) section,
this configuration should be held within a section called `custom_data` with a section-specific id.
The only difference is that no `data` subsection is given and a search pattern for the given id must
be supplied.

Search patterns are added [as with any other module](../getting_started/config.md#module-search-patterns).
Ensure that the search pattern key is the same as your `custom_data` section ID.

For example, a MultiQC config file could look as follows:

```yaml
# Other MultiQC config stuff here
custom_data:
  example_files:
    file_format: "tsv"
    section_name: "Coverage Decay"
    description: "This plot comes from files acommpanied by a multiqc_config.yaml file for configuration"
    plot_type: "linegraph"
    pconfig:
      id: "example_coverage_lineplot"
      title: "Coverage Decay"
      ylab: "X Coverage"
      ymax: 100
      ymin: 0
sp:
  example_files:
    fn: "example_files_*"
```

And work with the following data file:
`example_files_Sample_1.txt`:

```bash
0	98.22076066
1	97.96764159
2	97.78227175
3	97.61262195
# [...]
```

This kind of customisation should work with most Custom Content types.
For example, using an image called `some_science_mqc.jpeg` gives us a report section `some_science`,
which we can then add a nicer name and description to:

```yaml
custom_data:
  some_science:
    section_name: "Some real science"
    description: "This description comes from multiqc_config.yaml and helps to annotate the Custom Content image."
```

If no configuration is given, MultiQC will do its best to guess how to visualise
your data appropriately. To see examples of typical file structures which are understood, see the
[test data](https://github.com/MultiQC/test-data/tree/main/data/custom_content/no_config)
used to develop this code.

# Configuration

## Grouping sections and subsections

If you have multiple content types that you would like to group together with MultiQC sub-sections,
you can do so using the following keys:

```yaml
parent_id: custom_section
parent_name: "Some grouped data"
parent_description: "This parent section contains one or more sub-sections below it"
```

Any custom-content files that share the same `parent_id` will be grouped.

Note that some things, such as `parent_name` are taken from the first file that MultiQC finds
with this `parent_id`. So it's a good idea to specify this in every file.
`parent_description` and `extra` is taken from the first file where it is set.

:::warning
`parent_id` only works within Custom Content.
It is not currently possible to add custom content output into a report section
from a core MultiQC module.
:::

## Order of sections

If you have multiple different Custom Content sections, their order will be random
and may vary between runs. To avoid this, you can specify an order in your MultiQC
config as follows:

```yaml
custom_content:
  order:
    - first_cc_section
    - second_cc_section
```

Each section name should be the ID assigned to that section. You can explicitly set
this (see below), or the Custom Content module will automatically assign an ID.
To find out what your custom content section ID is, generate a report and click
the side navigation to your section. The browser URL should update and show something
that looks like this:

```
multiqc_report.html#my_cc_section
```

The section ID is the part after the `#` (`my_cc_section` in the above section).

Note that any Custom Content sections found that are _not_ specified in the config
will be placed at the top of the report.

## Section configuration

See below for how these config options can be specified (either within the data file
or in a MultiQC config file). All of these configuration parameters
are optional, and MultiQC will do its best to guess sensible defaults if they are
not specified.

All possible configuration keys and their default values are shown below:

```yaml
id: null # Unique ID for report section.
section_anchor: <id> # Used in report section #soft-links
section_name: <id> # Nice name used for the report section header
section_href: null # External URL for the data, to find more information
description: null # Introductory text to be printed under the section header
section_extra: null # Custom HTML to add after the section description
file_format: null # File format of the data (eg. csv / tsv)
plot_type:
  null # The plot type to visualise the data with.
  # generalstats | table | bargraph | linegraph | scatter | heatmap | beeswarm
pconfig: {} # Configuration for the plot.
```

:::info
Data types `generalstats` and `beeswarm` are _only_ possible by setting the above
configuration keys (these can't be guessed by data format).
:::

Note that any _custom content_ data found with the same section `id` will be merged
into the same report section / plot. The other section configuration keys are merged
for each file, with identical keys overwriting what was previously parsed.

This approach means that it's possible to have a single file containing data for multiple
samples, but it's also possible to have one file per sample and still have all of them
summarised.

:::note
If you're using `plot_type: 'generalstats'` then a report section will not be created and
most of the configuration keys above are ignored.
:::

## Plot configuration

Configuration of specific plots follows the same syntax as used when writing modules.
To find out more, please see the later docs. Specifically, the plot config docs for
[bar graphs](../development/plots.md#bar-graphs),
[line graphs](../development/plots.md#line-graphs),
[scatter plots](../development/plots.md#scatter-plots),
[tables](../development/plots.md#creating-a-table),
[beeswarm plots](../development/plots.md#beeswarm-plots-dot-plots) and
[heatmaps](../development/plots.md#heatmaps).

Wherever you see `pconfig`, any key can be used within the above syntax.

## Tricky extras

Because of the way this module works, there are a few specifics that can trip you up.
Most of these should probably be fixed one day. Feel free to ask for help on the [community forum](https://community.seqera.io/c/multiqc/6), or submit a pull request!
I'll try to keep a list here to help the wary...

### Differences between Tables and General Stats

Although they're both tables, note that general stats configures columns with a list
in the `pconfig` scope (see above example). Files that are just tables use `headers` instead.

### First columns in tables are special

The first column in every table is reserved for the sample name. As such, it shouldn't contain data.
All header configuration will be ignored for the first column. The only exception is name:
this can be tweaked using the somewhat tricky `col1_header` field in the `pconfig` scope (see table docs).
Alternatively, you can customise the column name by including a 'header row' in the first line of the `tsv`
or `csv` itself specifying the column names, with the first column with the name of your choice, and
subsequent columns including the key(s) defined in the header.

## Linting

MultiQC has been developed to be as forgiving as possible and will handle lots of
invalid or ignored configurations. This is useful for most users but can make life
difficult when getting MultiQC to work with a new custom content format.

To help with this, you can run MultiQC with the `--strict` flag, which will give
explicit warnings about anything that is not optimally configured. For example:

```bash
multiqc --strict test-data
```

You can alternatively enable the strict mode by setting the environment variable
`MULTIQC_STRICT`, or by setting it into the [config](http://multiqc.info/docs/#configuring-multiqc): `strict: true`.

# Examples

Probably the best way to get to grips with Custom Content is to see some examples.
The MultiQC automated testing runs with a
[number of different files](https://github.com/MultiQC/test-data/tree/main/data/custom_content)
which you can look through for inspiration.

For example, to see a file which generates a table in a report by itself, you can
have a look at `embedded_config/table_headers_mqc.txt` ([link](https://github.com/MultiQC/test-data/blob/main/data/custom_content/embedded_config/table_headers_mqc.txt)).
