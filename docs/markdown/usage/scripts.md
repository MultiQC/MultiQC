---
title: Using MultiQC in scripts
description: Importing MultiQC as a library in scripts and notebooks
---

# Using MultiQC within scripts

Even though the primary way to run MultiQC is as a command line, it can also be imported
like a Python module in order to build the report interactively,
such as in custom Python scripts or in a Jupyter notebook environment
(See an [example notebook](https://seqera.io/examples/jupyter/notebook.html)).

MultiQC provides a set of commands to iteratively parse logs and add sections to a report.
All of them are available via importing MultiQC as a module:

```python
import multiqc
```

## Parse logs

Find files that MultiQC recognizes in `analysis_dir` and parse them, without generating a report.
Data can be accessed with other methods: `list_modules`, `list_plots`, etc.

```python
def parse_logs(*analysis_dir, **kwargs)
```

Parameters:

- `analysis_dir`: Path(s) to search for files to parse
- `verbose`: Print more information to the console
- `file_list`: Supply a file containing a list of file paths to be searched, one per row
- `prepend_dirs`: Prepend directory to sample names
- `dirs_depth`: Prepend n directories to sample names. Negative number to take from start of path
- `fn_clean_sample_names`: Do not clean the sample names (leave as full file name)
- `require_logs`: Require all explicitly requested modules to have log files. If not, MultiQC will exit with an error
- `use_filename_as_sample_name`: Use the log filename as the sample name
- `strict`: Don't catch exceptions, run additional code checks to help development
- `quiet`: Only show log warnings
- `no_ansi`: Disable coloured log output
- `profile_runtime`: Add analysis of how long MultiQC takes to run to the report
- `no_version_check`: Disable checking the latest MultiQC version on the server
- `ignore`: Ignore analysis files
- `ignore_samples`: Ignore sample names
- `run_modules`: Use only this module. Can specify multiple times
- `exclude_modules`: Do not use this module. Can specify multiple times
- `config_files`: Specific config file to load, after those in MultiQC dir / home dir / working dir
- `module_order`: Names of modules in order of precedence to show in report
- `extra_fn_clean_exts`: Extra file extensions to clean from sample names
- `extra_fn_clean_trim`: Extra strings to clean from sample names
- `preserve_module_raw_data`: Preserve raw data from modules in the report - besides plots. Useful to use
  later interactively. Defaults to `True`. Set to `False` to save memory.

**Examples**

Parse logs found in the `data` directory.

```python
multiqc.parse_logs('data')
```

Parse logs found in the `data/fastp` directory, the `data/SAMPLE1.cutadapt.log` file,
and a `data_mqc.tsv` MultiQC [custom content](../custom_content/index.md) file.

```python
multiqc.parse_logs('data/fastp', 'data/SAMPLE1.cutadapt.log', "data_mqc.tsv")
```

Parse logs found in the `data` directory for only the specified modules, and use
and [additional pattern](../getting_started/config.md#cleaning-extensions) to clean sample names.

```python
multiqc.parse_logs(
    'data',
    run_modules=["fastp", "spades", "quast", "pangolin"],
    extra_fn_clean_exts=[".unclassified"],
)
```

Parse logs found in the `data` directory and run FastQC module twice for two sets of files - raw and trimmed reads - according to the provided path pattern (see [Order of modules](../reports/customisation.md#order-of-modules) for details).

```python
multiqc.parse_logs(
    'data',
    module_order=[
        dict(
            fastqc=dict(
                name="FastQC (trimmed)",
                anchor="fastqc_trimmed",
                path_filters=["*_1_trimmed_fastqc.zip"],
            )
        ),
        dict(
            quast=dict(
                name="FastQC (raw)",
                anchor="fastqc_raw",
                path_filters=["*_1_fastqc.zip"],
            )
        ),
    ],
)
```

MultiQC v1.29 and higher generates `BETA-multiqc.parquet` file in `multiqc_data` output directory. You can pass that file to `parse_logs`, and it will load that previous MultiQC run into memory.

Example:

```python
multiqc.parse_logs('multiqc_data/BETA-multiqc.parquet')
```

## List what's loaded

Return `list` of the modules that have been loaded, ordered according to config:

```python
def list_modules() ‑> list[str]
```

Return `dict` of plot names that have been loaded, indexed by module name and section:

```python
def list_plots() ‑> dict[str, list[str | dict[str, str]]]]
```

Example:

```python
multiqc.list_plots()
{'fastp': ['Filtered Reads',
  'Insert Sizes',
  {'Sequence Quality': ['Read 1: Before filtering',
    'Read 1: After filtering',
    'Read 2: Before filtering',
    'Read 2: After filtering']},
  {'GC Content': ['Read 1: Before filtering',
    'Read 1: After filtering',
    'Read 2: Before filtering',
    'Read 2: After filtering']},
  {'N content': ['Read 1: Before filtering',
    'Read 1: After filtering',
    'Read 2: Before filtering',
    'Read 2: After filtering']}]}
```

Return `list` of clean sample names that have loaded data:

```python
def list_samples() ‑> list[str]
```

Example:

```python
multiqc.list_samples()
['SAMPLE1_PE', 'SAMPLE2_PE']
```

Return `list` of found log files corresponding to the loaded data:

```python
def list_data_sources() ‑> list[str]
```

Example:

```python
multiqc.list_data_sources()
['data/SAMPLE1_PE.fastp.json', 'data/SAMPLE2_PE.fastp.json']
```

## Access loaded data

There are several methods to access the data loaded by `parse_logs`.

Return parsed module data, indexed (if available) by data key, then by sample. Module is either
the module name, or the anchor:

```python
def get_module_data(module: str = None, sample: str = None, key: str = None) ‑> dict
```

The function takes data from `report.saved_raw_data`, which populated by self.write_data_file() calls in individual modules.
This data is not necessarily normalized, e.g. numbers can be strings or numbers, depending on the individual module behaviour.

Example:

```python
> multiqc.get_module_data(module="fastp", sample="SAMPLE1_PE")["summary"]
{'fastp_version': '0.23.2',
 'sequencing': 'paired end (301 cycles + 301 cycles)',
 'before_filtering': {'total_reads': 55442,
  'total_bases': 16571632,
  'q20_bases': 16267224,
  'q30_bases': 15853021,
  'gc_content': 0.38526},
 'after_filtering': {'total_reads': 48270,
  'total_bases': 14363465,
  'q20_bases': 14323363,
  'q30_bases': 14199841,
  'gc_content': 0.383991}}
```

Similarly, return parsed general stats data, indexed by sample, then by data key. If sample is specified, return only data for that sample.

```python
def get_general_stats_data(sample: str = None) ‑> dict
```

## Adding custom content

You can also custom section to the report by subclassing from `multiqc.BaseMultiqcModule`. This can be used to add a custom table or other content.

**Example**

Create a table (see [plotting](../development/plots.md#creating-a-table) for more detail) and add it to the report.

```python
import multiqc
from multiqc.plots import table

plot = table.plot(
    data=...,
    headers=...,
    pconfig={
        "id": "my_metrics_table",
        "title": "My metrics",
    },
)
module = multiqc.BaseMultiqcModule(
    name="my-module",
    anchor="custom_data",
)
module.add_section(
    plot=plot,
    name="My metrics",
    anchor="my_metrics_section",
    description=...,
)
multiqc.report.modules.append(module)
```

## Get plot object

Get a plot object for a specific module and section. For list of available plots, use `multiqc.list_plots`.

```python
def get_plot(module: str, section: str) -> Plot
```

**Examples**

Get plot object for the "GC Content" plot in the "fastp" module.

```python
plot = multiqc.get_plot("fastp", "GC Content")
```

Get plot object for the "Number of Contigs" plot in the "QUAST" module.

```python
plot = multiqc.get_plot("QUAST", "Number of Contigs")
```

## Show plot

Prepare plot to be shown in the notebook cell.

```python
class Plot:
    def show(self, dataset_id: int | str = 0, flat=False, **kwargs)
```

Parameters:

- `dataset_id`: Dataset label, in case if plot has several tabs
- `flat`: Show plot as static images without any interactivity
- `kwargs`: Additional arguments passed to the plot

**Examples**

Create a bar graph and show it in the notebook cell:

```python
from multiqc.plots import bargraph
plot = bargraph.plot(...)
display(plot.show(violin=True))
```

Get "fastp GC Content" plot and show it in the notebook cell. Since it has multiple
tabs, we can select which tab to show with the `dataset_id` option (defaults to the first tab):

```python
plot = multiqc.get_plot("fastp", "GC Content")
display(plot.show(dataset_id="Read 2: Before filtering"))
```

Shows Samtools alignment stats as a violin plot. Use flat image without interactivity.

```python
plot = multiqc.get_plot("Flagstat", "Alignment stats")
display(plot.show("Read counts", violin=True, flat=True))
```

:::note
Calling the notebook's built-in `display` function is optional when the `show` call is the last line in your cell.
:::

## Save plot to file

Similarly, you can save plot to a file instead of showing it in a notebook.

```python
class Plot:
    def save(self, filename, dataset_id: int | str = 0, **kwargs)
```

Parameters:

- `filename`: Path to save the plot
- `dataset_id`: Dataset label, in case if plot has several tabs
- `kwargs`: Additional arguments passed to the plot

**Examples**:

Save the "Number of Contigs" plot for the QUAST module to a file.

```python
plot = multiqc.get_plot("QUAST", "Number of Contigs")
plot.save("quast_contigs.html")
```

Save the GC Content plot for the dataset labeled "Read 2: Before filtering" to a file,
make it flat.

```python
plot = multiqc.get_plot("fastp", "GC Content")
plot.save(
    "fastp_gc_content.png",
    dataset_id="Read 2: Before filtering",
)
```

Save Samtools alignment stats as a violin plot to a file.

```python
plot = multiqc.get_plot("Flagstat", "Alignment stats")
plot.save(
    "flagstat_alignment_stats.html",
    dataset="Read counts",
    violin=True,
)
```

## Writing report

Render HTML from parsed module data, and write a report along with auxiliary data files to disk.

```python
def write_report(**kwargs)
```

Parameters:

- `title`: Report title. Printed as page header, used for filename if not otherwise specified
- `report_comment`: Custom comment, will be printed at the top of the report
- `template`: Report template to use
- `output_dir`: Create report in the specified output directory
- `filename`: Report filename. Use 'stdout' to print to standard out
- `make_data_dir`: Force the parsed data directory to be created
- `data_format`: Output parsed data in a different format
- `zip_data_dir`: Compress the data directory
- `force`: Overwrite existing report and data directory
- `make_report`: Generate the report HTML. Defaults to `True`, set to `False` to only export data and plots
- `export_plots`: Export plots as static images in addition to the report
- `plots_force_flat`: Use only flat plots (static images)
- `plots_force_interactive`: Use only interactive plots (in-browser Javascript)
- `strict`: Don't catch exceptions, run additional code checks to help development
- `development`: Development mode. Do not compress and minimise JS, export uncompressed plot data
- `make_pdf`: Create PDF report. Requires Pandoc to be installed
- `no_megaqc_upload`: Don't upload generated report to MegaQC, even if MegaQC options are found
- `quiet`: Only show log warnings
- `verbose`: Print more information to the console
- `no_ansi`: Disable coloured log output
- `profile_runtime`: Add analysis of how long MultiQC takes to run to the report
- `no_version_check`: Disable checking the latest MultiQC version on the server
- `run_modules`: Use only these modules
- `exclude_modules`: Do not use these modules
- `config_files`: Specific config file to load, after those in MultiQC dir / home dir / working dir
- `custom_css_files`: Custom CSS files to include in the report
- `module_order`: Names of modules in order of precedence to show in report

Example:

```python
multiqc.write_report(
    force=True,
    output_dir="my_multiqc_report",
    title="My Report",
    filename="report.html",
)
```

## Load config from file

Load config on top of the current config from a MultiQC config file.

```python
def load_config(config_file: str | Path)
```

## Reset session

Reset the report to start fresh. Drops all previously parsed and loaded data.

```python
def reset()
```
