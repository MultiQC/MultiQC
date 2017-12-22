# Customising Reports
MultiQC offers a few ways to customise reports to easily add your own
branding and some additional report-level information. These features
are primarily designed for core genomics facilities.

Note that much more extensive customisation of reports is possible using
[custom templates](http://multiqc.info/docs/#writing-new-templates).

## Titles and introductory text
You can specify a custom title for the report using the `-i`/`--title`
command line option. The `-b`/`--comment` option can be used to add a
longer comment to the top of the report at run time.

You can also specify the title and comment, as well as
a subtitle and the introductory text in your config file:

```yaml
title: "My Title"
subtitle: "A subtitle to go underneath in grey"
intro_text: "MultiQC reports summarise analysis results."
report_comment: "This is a comment about this report."
```

Note that if `intro_text` is `None` the template will display the default
introduction sentence. Set this to `False` to hide this, or set it to a
string to use your own text.

## Report Logo
To add your own custom logo to reports, you can add the following
three lines to your MultiQC configuration file:

```yaml
custom_logo: '/abs/path/to/logo.png'
custom_logo_url: 'https://www.example.com'
custom_logo_title: 'Our Institute Name'
```
Only `custom_logo` is needed. The URL will make the logo open up
a new web browser tab with your address and the title sets the mouse
hover title text.

## Project level information
You can add custom information at the top of reports by adding key:value
pairs to the config option `report_header_info`. Note that if you have
a file called `multiqc_config.yaml` in the working directory, this will
automatically be parsed and added to the config. For example, if you have
the following saved:

```yaml
report_header_info:
    - Contact E-mail: 'phil.ewels@scilifelab.se'
    - Application Type: 'RNA-seq'
    - Project Type: 'Application'
    - Sequencing Platform: 'HiSeq 2500 High Output V4'
    - Sequencing Setup: '2x125'
```

Then this will be displayed at the top of reports:

![report project info](images/report_proj_info.png)

Note that you can also specify a path to a config file using `-c`.

## Bulk sample renaming
Although it is possible to rename samples manually and in bulk using the
[report toolbox](#renaming-samples), it's often desirable to embed such renaming patterns
into the report so that they can be shared with others. For example, a typical case could be
for a sequencing centre that has internal sample IDs and also user-supplied sample names.
Or public sample identifiers such as SRA numbers as well as more meaningful names.

It's possible to supply a file with one or more sets of sample names using the `--sample-names`
command line option. This file should be a tab-delimited file with a header row (used for
the report button labels) and then any number of renamed sample identifiers. For example:

```
MultiQC Names	Proper Names	AWESOME NAMES
SRR1067503_1	Sample_1	MYBESTSAMP_1
SRR1067505_1	Sample_2	MYBESTSAMP_2
SRR1067510_1	Sample_3	MYBESTSAMP_3
```

If supplied, buttons will be generated at the top of the report with your labels.
Clicking these will populate and apply the Toolbox renaming panel.

> **NB:** Sample renaming works with partial substrings - these will be replaced!

It's also possible to supply such renaming patterns within a config file (useful if you're
already generating a config file for a run). In this case, you need to set the variables
`sample_names_rename_buttons` and `sample_names_rename`. For example:

```yaml
sample_names_rename_buttons:
    - "MultiQC Names"
    - "Proper Names"
    - "AWESOME NAMES"
sample_names_rename:
    - ["SRR1067503_1", "Sample_1", "MYBESTSAMP_1"]
    - ["SRR1067505_1", "Sample_2", "MYBESTSAMP_2"]
    - ["SRR1067510_1", "Sample_3", "MYBESTSAMP_3"]
```

## Module and section comments
Sometimes you may want to add a custom comment above specific sections in the report. You can
do this with the config option `section_comments` as follows:

```yaml
section_comments:
    featurecounts: 'This comment is for a module header, but should still work'
    star_alignments: 'This new way of commenting above sections is **awesome**!'
```

Comments can be written in Markdown. The `section_comments` keys should correspond to the HTML IDs
of the report section. You can find these by clicking on a navigation link in the report and seeing
the `#section_id` at the end of the browser URL.

## Order of modules
By default, modules are included in the report as in the order specified in `config.module_order`.
Any modules found which aren't in this list are appended at the top of the report.

#### Top modules
To specify certain modules that should always come at the top of the report, you can configure
`config.top_modules` in your MultiQC configuration file. For example, to always have the FastQC
module at the top of reports, add the following to your `~/.multiqc_config.yaml` file:

```yaml
top_modules:
    - 'fastqc'
```

#### Running modules multiple times
A module can be specified multiple times in either `config.module_order` or `config.top_modules`,
causing it to be run multiple times. By itself you'll just get two identical report sections.
However, you can also supply configuration options to the modules as follows:

```yaml
top_modules:
    - moduleName:
        name: 'Module (filtered)'
        info: 'This section shows the module with different files'
        path_filters:
            - '*_special.txt'
            - '*_others.txt'
    - moduleName:
        name: 'Module (all)'
```
These overwrite the defaults that are hardcoded in the module code. `path_filters` is the
exception, which filters the file searches for a given list of glob filename patterns.
The available options are:

* `name`: Section name
* `anchor`: Section report ID
* `target`: Intro link text
* `href`: Intro link URL
* `info`: Intro text
* `extra`: Additional HTML after intro.
* `custom_config`: Custom module-level settings. Translated into `config.moduleName`, but specifically for this section.

For example, to run the FastQC module twice, before and after adapter trimming, you could
use the following config:

```yaml
module_order:
    - fastqc:
        name: 'FastQC (trimmed)'
        info: 'This section of the report shows FastQC results after adapter trimming.'
        target: ''
        path_filters:
            - '*_1_trimmed_fastqc.zip'
    - cutadapt
    - fastqc:
        name: 'FastQC (raw)'
        path_filters:
            - '*_1_fastqc.zip'
```

Note that if you change the `name` then you will get multiples of columns in the
_General Statistics_ table. If unchanged, the topmost module may overwrite output from
the first iteration.

> NB: Currently, you can not list a module name in both `top_modules` and `module_order`.
> Let me know if this is a problem..


## Order of sections
Sometimes it's desirable to customise the order of specific sections in a report, independent of
module execution. For example, the `custom_content` module can generate multiple sections from
different input files.

To do this, follow a link in a report navigation to skip to the section you want to move (must
be a major section header, not a subheading). Find the ID of that section by looking at the URL.
For example, clicking on _FastQC_ changes the URL to `multiqc_report.html#fastqc` -  the ID is
the text after (not including) the `#` symbol.

Next, specify the `report_section_order` option in your MultiQC config file. Section in
the report are given a number ranging from 10 (section at bottom of report), incrementing by +10
for each section. You can change this number (eg. a very low number to always get at the bottom
of the report or very high to always be at the top), or you can move a section to before or after
another existing section (has no effect if the other named ID is not in the report).

For example, add the following to your MultiQC config file:

```yaml
report_section_order:
    section1:
        order: -1000
    section2:
        before: 'othersection'
    section3:
        after: 'diffsection'
```

## Customising tables

### Hiding columns
Report tables such as the General Statistics table can get quite wide. To help with this,
columns in the report can be hidden. Some MultiQC modules include columns which are hidden
by default, others may be uninteresting to some users.

To allow customisation of this behaviour, the defaults can be changed by adding to your
MultiQC config file. This is done with the `table_columns_visible` value. Open a MultiQC
report and click _Configure Columns_ above a table. Make a note of the _Group_ and _ID_
for the column that you'd like to alter. For example, to make the `% Duplicate Reads`
column from FastQC hidden by default, the _Group_ is `FastQC` and the _ID_ is
`percent_duplicates`. These are then added to the config as follows:

```yaml
table_columns_visible:
    FastQC:
        percent_duplicates: False
```

Note that you can set these to `True` to show columns that would otherwise be hidden
by default.

### Column order
In the same way, you can force a column to appear at the start or end of the table, or
indeed impose a custom ordering on all the columns, by setting the `table_columns_placement`.
High values push columns to the right hand side of the table and low to the left. The default
value is 1000. For example:

```yaml
table_columns_placement:
    Samtools:
        reads_mapped: 900
        properly_paired: 1010
        secondary: 1020
```

In this case, since the default placement weighting is `1000`, the `reads_mapped` will end up as the
leftmost column and the other two will and up as the final columns on the right of the table.

### Conditional formatting
It's possible to highlight values in tables based on their value. This is done using the `table_cond_formatting_rules` config setting. Rules can be applied to every table column, or to specific columns only, using that column's unique ID.

The default rules are as follows:

```yaml
table_cond_formatting_rules:
    all_columns:
        pass:
            - s_eq: 'pass'
            - s_eq: 'true'
        warn:
            - s_eq: 'warn'
            - s_eq: 'unknown'
        fail:
            - s_eq: 'fail'
            - s_eq: 'false'
```

These make any table cells that match the string `pass` or `true` have text with a green background, orange for `warn`, red for `fail` and so on. There can be multiple tests for each style of formatting - if there is a match for any, it will be applied. The following comparison operators are available:

* `s_eq` - String exactly equals (case insensitive)
* `s_contains` - String contains (case insensitive)
* `s_ne` - String does not equal (case insensitive)
* `eq` - Value equals
* `ne` - Value does not equal
* `gt` - Value is greater than
* `lt` - Value is less than

To have matches for a specific column, use that column's ID instead of `all_columns`. For example:

```yaml
table_cond_formatting_rules:
    mqc-generalstats-uniquely_mapped_percent:
        pass:
            - gt: 80
        warn:
            - lt: 80
        fail:
            - lt: 70
```

Note that the formatting is done in a specific order - `pass`/`warn`/`fail` by default, so that anything matching both `warn` and `fail` will be formatted as `fail` for example. This can be customised with `table_cond_formatting_colours` (see below).

To find the unique ID for your column, right click a table cell in a report and inspect it's HTML (_Inpsect_ in Chrome). It should look something like `<td class="data-coloured mqc-generalstats-Assigned">`, where the `mqc-generalstats-Assigned` bit is the unique ID.

> I know this isn't the same method of IDs as above and isn't super easy to do. Sorry!

It's possible to highlight matches in any number of colours. MultiQC comes with the following defaults:

```yaml
table_cond_formatting_colours:
    - blue: '#337ab7'
    - lbue: '#5bc0de'
    - pass: '#5cb85c'
    - warn: '#f0ad4e'
    - fail: '#d9534f'
```

These can be overridden or added to with any string / CSS hex colour combinations you like. You can generate hex colour codes with lots of tools, for example http://htmlcolorcodes.com/

Note that the different sets of rules are formatted in order. So if a value matches both `pass` and `fail` then it will be formatted as a `fail`


## Number base (multiplier)
To make numbers in the General Statistics table easier to read and compare quickly,
MultiQC sometimes divides them by one million (typically read counts). If your
samples have very low read counts then this can result in the table showing
counts of `0.0`, which isn't very helpful.

To change this behaviour, you can customise three config variables in your MultiQC
config. The defaults are as follows:
```yaml
read_count_multiplier: 0.000001
read_count_prefix: 'M'
read_count_desc: 'millions'
```

So, to show thousands of reads instead of millions, change these to:
```yaml
read_count_multiplier: 0.001
read_count_prefix: 'K'
read_count_desc: 'thousands'
```

The same options are also available for numbers of base pairs:
```yaml
base_count_multiplier: 0.000001
base_count_prefix: 'Mb'
base_count_desc: 'millions'
```


## Number formatting
By default, the interactive HighCharts plots in MultiQC reports use spaces for thousand
separators and points for decimal places (_e.g._ `1 234 567.89`). Different countries
have different preferences for this, so you can customise the two using a couple of
configuration parameters - `decimalPoint_format` and `thousandsSep_format`.

For example, the following config would result in the following alternative
number formatting: `1234567,89`.
```yaml
decimalPoint_format: ','
thousandsSep_format: ''
```

This formatting currently only applies to the interactive charts. It may be extended
to apply elsewhere in the future (submit a new issue if you spot somewhere where you'd like it).


## Troubleshooting
One tricky bit that caught me out whilst writing this is the different type casting
between Python, YAML and Jinja2 templates. This is especially true when using an
empty variable:
```python
# Python
my_var = None
```
```yaml
# YAML
my_var: null
```
```python
# Jinja2
if myvar is none # Note - Lower case!
```
