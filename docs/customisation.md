# Customising Reports
MultiQC offers a few ways to customise reports to easily add your own
branding and some additional report-level information. These features
are primarily designed for core genomics facilities.

Note that much more extensive customisation of reports is possible using
[custom templates](#writing-new-templates).

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
    Contact E-mail: 'phil.ewels@scilifelab.se'
    Application Type: 'RNA-seq'
    Project Type: 'Application'
    Sequencing Platform: 'HiSeq 2500 High Output V4'
    Sequencing Setup: '2x125'
```

Then this will be displayed at the top of reports:

![report project info](images/report_proj_info.png)

Note that you can also specify a path to a config file using `-c`.


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

