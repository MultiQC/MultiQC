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

## Order of modules
By default, modules are included in the report as in the order specified in `config.module_order`.
Any modules found which aren't in this list are appended at the top of the report. To specify
certain modules that should always come at the top of the report, you can configure `config.top_modules`
in your MultiQC configuration file. For example, to always have the FastQC module at the top
of reports, add the following to your `~/.multiqc_config.yaml` file:

```yaml
top_modules:
    - 'fastqc'
```

## Customising tables
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
