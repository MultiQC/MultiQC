
> **WARNING** - This feature is new and is very much in a beta status.
> It is expected to be further developed in future releases, which may break backwards
> compatibility. There are also probably quite a few bugs. Use at your own risk!
> Please report bugs or missing functionality as a
> [new GitHub issue](https://github.com/ewels/MultiQC/issues/new).

# Introduction
Bioinformatics projects often include non-standardised analyses, with results from custom
scripts or in-house packages. It can be frustrating to have a MultiQC report describing
results from 90% of your pipeline but missing the final key plot. To help with this,
MultiQC has a special _"custom content"_ module.

Custom content parsing is a little more restricted than standard modules. Specifically:

* Only one plot per section is possible
* Plot customisation is more limited

All plot types can be generated using custom content - see the
[test files](https://github.com/ewels/MultiQC_TestData/tree/master/data/custom_content)
for examples of how data should be structured.

# Configuration
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
id: null                # Unique ID for report section.
section_anchor: <id>    # Used in report section #soft-links
section_name: <id>      # Nice name used for the report section header
section_href: null      # External URL for the data, to find more information
description: null       # Introductory text to be printed under the section header
file_format: null       # File format of the data (typically csv / tsv - see below for more information)
plot_type: null         # The plot type to visualise the data with.
                        # - Possible options: generalstats | table | bargraph | linegraph | scatter | heatmap | beeswarm
pconfig: {}             # Configuration for the plot. See http://multiqc.info/docs/#plotting-functions
```

Note that any _custom content_ data found with the same section `id` will be merged
into the same report section / plot. The other section configuration keys are merged
for each file, with identical keys overwriting what was previously parsed.

This approach means that it's possible to have a single file containing data for multiple
samples, but it's also possible to have one file per sample and still have all of them
summarised.

If you're using `plot_type: 'generalstats'` then a report section will not be created and
most of the configuration keys above are ignored.

Data types `generalstats` and `beeswarm` are _only_ possible by setting the above
configuration keys (these can't be guessed by data format).

# Data formats
MultiQC can parse custom data from a few different sources, in a number of different
formats. Which one you use depends on how the data is being produced.

A quick summary of which approach to use looks something like this:

* Additional data when already using custom MultiQC config files
  * _Data as part of MultiQC config_
* Data specifically for MultiQC from a custom script
  * _MultiQC-specific data file_
* Data from a custom script which is also used by other processes
  * _Separate configuration and data files_
  * Add `_mqc.txt` to filename and hope that MultiQC guesses correctly
* Anything more complicated, or data from a released tool
  * Write a proper MultiQC module instead.

For more complete examples of the data formats understood by MultiQC, please see the
[`data/custom_content`](https://github.com/ewels/MultiQC_TestData/tree/master/data/custom_content)
directory in the [MultiQC_TestData](https://github.com/ewels/MultiQC_TestData)
GitHub repository.

## Data from a released tool
If your data comes from a released bioinformatics tool, you shouldn't be using this
feature of MultiQC! Sure, you can probably get it to work, but it's better if a
fully-fledged core MultiQC module is written instead. That way, other users of MultiQC
can also benefit from results parsing.

Note that proper MultiQC modules are more robust and powerful than this custom-content
feature. You can also [write modules](http://multiqc.info/docs/#writing-new-modules)
in [MultiQC plugins](http://multiqc.info/docs/#multiqc-plugins) if they're not suitable for
general release.

## Data as part of MultiQC config
If you are already using a MultiQC config file to add data to your report (for example,
[titles / introductory text](http://multiqc.info/docs/#customising-reports)), you can
give data within this file too. This can be in any MultiQC config file (for example,
passed on the command line with `-c my_yaml_file.yaml`). This is useful as you can
keep everything contained within a single file (including stuff unrelated to this
specific _custom content_ feature of MultiQC).

If you're not using this file for other MultiQC configuration, you're probably
better off using a stand-alone YAML file (see [section below](#multiqc-specific-data-file)).

To be understood by MultiQC, the `custom_data` key must be found. This must contain
a section with a unique id, specific to your new report section. This in turn
must contain a section called `data`. Other configuration keys can be held alongside this.
For example:

```yaml
# Other MultiQC config stuff here
custom_data:
    my_data_type:
        id: 'mqc_config_file_section'
        section_name: 'My Custom Section'
        description: 'This data comes from a single multiqc_config.yaml file'
        plot_type: 'bargraph'
        pconfig:
            id: 'barplot_config_only'
            title: 'MultiQC Config Data Plot'
            ylab: 'Number of things'
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
        plot_type: 'generalstats'
        pconfig:
            - col_1:
                max: 100
                min: 0
                scale: 'RdYlGn'
                suffix: '%'
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
> **Note:** Use a **list** of headers in `pconfig` (keys prepended with `-`) to specify the order
> of columns in the table.

See the [general statistics docs](http://multiqc.info/docs/#step-4-adding-to-the-general-statistics-table)
for more information about configuring data for the General Statistics table.

## MultiQC-specific data file
If you can choose exactly how your data output looks, then the easiest way to parse it
is to use a MultiQC-specific format. If the filename ends in `*_mqc.(yaml|json|txt|csv|out)`
then it will be found by any standard MultiQC installation with no additional customisation
required (v0.9 onwards).

These files contain configuration information specifying how the data should be parsed, along
side the data. If using YAML, this looks just the same as if in a MultiQC config file (see above),
but without having to be within a `custom_data` section:

```yaml
id: 'my_pca_section'
section_name: 'PCA Analysis'
description: 'This plot shows the first two components from a principal component analysis.'
plot_type: 'scatter'
pconfig:
    id: 'pca_scatter_plot'
    title: 'PCA Plot'
    xlab: 'PC1'
    ylab: 'PC2'
data:
    sample_1: {x: 12, y: 14}
    sample_2: {x: 8, y: 6 }
    sample_3: {x: 5, y: 11}
    sample_4: {x: 9, y: 12}
```

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

If you want the data to be easy to use with other tools, you can also use
comma-separated or tab-separated file. To customise plot output, include commented
header lines with plot configuration in YAML format:

```bash
# title: 'Output from my script'
# description: 'This output is described in the file header. Any MultiQC installation will understand it without prior configuration.'
# section: 'Custom Data File'
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

If no configuration is given, MultiQC will do its best to guess how to visualise your data appropriately.
To see examples of typical file structures which are understood, see the
[test data](https://github.com/ewels/MultiQC_TestData/tree/master/data/custom_content/no_config)
used to develop this code. Something will be probably be shown, but it may produce unexpected results.

## Separate configuration and data files
It's not always possible or desirable to include MultiQC configuration within a data file.
If this is the case, you can add to the MultiQC configuration to specify how input files
should be parsed.

As described in the above [_Data as part of MultiQC config_](#data-as-part-of-multiqc-config) section,
this configuration should be held within a section called `custom_data` with a section-specific id.
The only difference is that no `data` subsection is given and a search pattern for the given id must
be supplied.

Search patterns are added [as with any other module](http://multiqc.info/docs/#module-search-patterns).
Ensure that the search pattern key is the same as your `custom_data` section ID.

For example:
```yaml
# Other MultiQC config stuff here
custom_data:
    example_files:
        file_format: 'tsv'
        section_name: 'Coverage Decay'
        description: 'This plot comes from files acommpanied by a mutliqc_config.yaml file for configuration'
        plot_type: 'linegraph'
        pconfig:
            id: 'example_coverage_lineplot'
            title: 'Coverage Decay'
            ylab: 'X Coverage'
            ymax: 100
            ymin: 0
sp:
    example_files:
        fn: 'example_files_*'
```
A data file within the MultiQC search directories could then simply look like this:

`example_files_Sample_1.txt`:
```
0	98.22076066
1	97.96764159
2	97.78227175
3	97.61262195
[...]
```

As mentioned above - if no configuration is given, MultiQC will do its best to guess how to visualise
your data appropriately. To see examples of typical file structures which are understood, see the
[test data](https://github.com/ewels/MultiQC_TestData/tree/master/data/custom_content/no_config)
used to develop this code.

