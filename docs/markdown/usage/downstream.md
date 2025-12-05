---
title: Downstream analysis
description: How to use MultiQC raw data outputs
---

# Downstream analysis

Whilst MultiQC is typically used as a final reporting step in an analysis, it can also be used as an intermediate in your analysis.

MultiQC saves a directory of machine-readable outputs called `multiqc_data/`. In here there are files from each module and table, as well as a verbose `multiqc.log` file and, a `multiqc.parquet` file that contains all the intermediate data and metadata needed to regenereate a report.

Most of these files are tab-separated `.tsv` files by default, but you can choose to have them as JSON, YAML if you prefer with the `-k`/`--data-format` flag or the `data_format` option in a config file.

These files can be useful as MultiQC essentially standardises the outputs from a lot of different tools.
Typical usage of MultiQC outputs could be filtering of large datasets (eg. single-cell analysis) or trend-monitoring of repeated runs.

## Programmatic HTML Report Access

MultiQC provides clean programmatic access to HTML report content, making it easy to integrate MultiQC into web applications, Jupyter notebooks, Streamlit dashboards, and other interactive environments.

### Method 1: Using `multiqc.run()` with `return_html=True`

```python
import multiqc
from multiqc.core.update_config import ClConfig

result = multiqc.run(
    "./analysis_data",
    cfg=ClConfig(quiet=True, no_data_dir=True),
    return_html=True
)

if result.sys_exit_code == 0 and result.html_content:
    html_content = result.html_content
    print(f"Generated {len(html_content)} characters of HTML")
    # Use html_content in your application
else:
    print(f"MultiQC failed: {result.message}")
```

### Method 2: Using `multiqc.write_report()` with `return_html=True`

```python
import multiqc

# Parse the analysis files first
multiqc.parse_logs("./analysis_data", quiet=True)

# Generate HTML report
html_content = multiqc.write_report(
    title="My Quality Control Report",
    quiet=True,
    return_html=True
)

if html_content:
    # Use html_content in Streamlit, Jupyter, Flask, etc.
    pass
```

### Integration Examples

#### Streamlit Dashboard

```python
import streamlit as st
import multiqc
from multiqc.core.update_config import ClConfig

st.title("Quality Control Dashboard")
uploaded_files = st.file_uploader("Upload analysis files", accept_multiple_files=True)

if uploaded_files:
    # Process files and generate report
    with tempfile.TemporaryDirectory() as temp_dir:
        # Save uploaded files to temp directory
        for file in uploaded_files:
            with open(os.path.join(temp_dir, file.name), "wb") as f:
                f.write(file.read())

        # Generate MultiQC report
        result = multiqc.run(temp_dir, cfg=ClConfig(quiet=True), return_html=True)

        if result.html_content:
            st.components.v1.html(result.html_content, height=800, scrolling=True)
```

#### Jupyter Notebook

```python
from IPython.display import HTML
import multiqc

result = multiqc.run("./experiment_data", return_html=True)
if result.html_content:
    display(HTML(result.html_content))
```

#### Flask Web Application

```python
from flask import Flask
import multiqc

app = Flask(__name__)

@app.route('/qc-report/<path:analysis_path>')
def generate_qc_report(analysis_path):
    result = multiqc.run(analysis_path, return_html=True)
    return result.html_content if result.html_content else "Error generating report"
```

### Notes

- The HTML content is fully self-contained with all CSS and JavaScript
- Use `quiet=True` and `no_data_dir=True` for cleaner programmatic usage
- Always check `result.sys_exit_code` and `result.html_content` before using the HTML
- The returned HTML can be several MB in size for large datasets

Below are a few tools that are specifically designed to work with MultiQC.
They are not created by or endorsed by the MultiQC author but may be helpful for your research.

## TidyMultiqc

- Homepage: [https://CRAN.R-project.org/package=TidyMultiqc](https://CRAN.R-project.org/package=TidyMultiqc)
- Source: [https://github.com/TMiguelT/TidyMultiqc](https://github.com/TMiguelT/TidyMultiqc)

Provides the means to convert `multiqc_data.json` files into `tidy` data frames for downstream analysis in R.

This analysis might involve cohort analysis, quality control visualisation, change-point detection, statistical process control, clustering, or any other type of quality analysis.

## Programmatic HTML Report Capture

MultiQC supports capturing the HTML report content programmatically, which is useful for integrating MultiQC into web applications, Jupyter notebooks, or other interactive environments like Streamlit dashboards.

### Method 1: Using filename='stdout' with Python API

```python
import sys
from io import StringIO
import multiqc
from multiqc.core.update_config import ClConfig

def get_multiqc_html(analysis_dir):
    """Capture MultiQC HTML report as a string variable"""
    # Save the original stdout
    original_stdout = sys.stdout
    captured_output = StringIO()

    try:
        # Redirect stdout to capture the HTML
        sys.stdout = captured_output

        # Configure MultiQC to output HTML to stdout
        cfg = ClConfig(
            filename="stdout",
            quiet=True,          # Suppress log messages
            no_data_dir=True     # Don't create data directory
        )

        # Run MultiQC
        result = multiqc.run(analysis_dir, cfg=cfg)

        if result.sys_exit_code == 0:
            return captured_output.getvalue()
        return None

    finally:
        # Always restore stdout
        sys.stdout = original_stdout
        captured_output.close()

# Usage example
html_content = get_multiqc_html("./analysis_data")
if html_content:
    print(f"Captured {len(html_content)} characters of HTML content")
    # Use html_content in your application
```

### Method 2: Using CLI with subprocess

```python
import subprocess
import sys

def get_multiqc_html_cli(analysis_dir):
    """Capture MultiQC HTML using CLI interface"""
    result = subprocess.run([
        sys.executable, "-m", "multiqc",
        analysis_dir,
        "--filename", "stdout",
        "--quiet",
        "--no-data-dir"
    ], capture_output=True, text=True)

    if result.returncode == 0:
        return result.stdout
    return None

# Usage
html_content = get_multiqc_html_cli("./data")
```

### Integration Examples

#### Streamlit Dashboard

```python
import streamlit as st

def generate_multiqc_report(data_dir):
    # Use Method 1 from above
    return get_multiqc_html(data_dir)

st.title("Quality Control Dashboard")

# File uploader
uploaded_files = st.file_uploader(
    "Upload analysis files",
    accept_multiple_files=True
)

if uploaded_files:
    # Process uploaded files and generate report
    html_report = generate_multiqc_report("./temp_data")

    if html_report:
        # Display the HTML report directly in Streamlit
        st.components.v1.html(html_report, height=800, scrolling=True)
    else:
        st.error("Failed to generate MultiQC report")
```

#### Jupyter Notebook

```python
from IPython.display import HTML

# Generate and display MultiQC report inline
html_content = get_multiqc_html("./experiment_data")
if html_content:
    display(HTML(html_content))
```

#### Flask Web Application

```python
from flask import Flask, render_template_string

app = Flask(__name__)

@app.route('/qc/<path:analysis_path>')
def quality_control_report(analysis_path):
    """Serve MultiQC report dynamically"""
    html_content = get_multiqc_html(analysis_path)

    if html_content:
        return html_content  # Return HTML directly
    else:
        return "Error generating report", 500

@app.route('/dashboard/<path:analysis_path>')
def embedded_dashboard(analysis_path):
    """Embed MultiQC report in a larger dashboard"""
    html_content = get_multiqc_html(analysis_path)

    dashboard_template = """
    <html>
    <head><title>Analysis Dashboard</title></head>
    <body>
        <div class="header">
            <h1>Lab Quality Control Dashboard</h1>
            <nav><!-- Navigation links --></nav>
        </div>
        <div class="multiqc-container">
            {{ multiqc_html|safe }}
        </div>
        <div class="footer">
            <p>Generated with MultiQC</p>
        </div>
    </body>
    </html>
    """

    return render_template_string(
        dashboard_template,
        multiqc_html=html_content
    )
```

### Notes

- The HTML content includes all CSS and JavaScript needed for interactive plots
- Use `quiet=True` to suppress log messages when capturing HTML
- Consider using `no_data_dir=True` if you don't need the raw data files
- The captured HTML is fully self-contained and can be saved, served, or embedded as needed
- This approach works with all MultiQC modules and configurations

## MegaQC

- Homepage: [https://megaqc.info](https://megaqc.info)
- Source: [https://github.com/ewels/MegaQC](https://github.com/ewels/MegaQC)

Started off by MultiQC author [@ewels](https://github.com/ewels/) this project has had further development by a team of several contributors. It is functional but still has several parts of its codebase that have never quite been finished.

MegaQC imports data from multiple MultiQC runs and provides an interface to explore this with an interactive web server using a database backend.
It can plot data over time, across runs and even has an interactive dashboard builder.
It's useful for anyone who wants to monitor MultiQC statistics (eg. clinical labs) or work interactively with large datasets (eg. single cell analysis).

## ChronQC

- Docs: [https://chronqc.readthedocs.io](https://chronqc.readthedocs.io)
- Source: [https://github.com/nilesh-tawari/ChronQC](https://github.com/nilesh-tawari/ChronQC)

ChronQC is a quality control (QC) tracking system for clinical implementation of next-generation sequencing (NGS). ChronQC generates time series plots for various QC metrics, which allows comparison of the current run to historical runs. ChronQC has multiple features for tracking QC data including Westgard rules for clinical validity, laboratory-defined thresholds, and historical observations within a specified period. Users can record their notes and corrective actions directly onto the plots for long-term recordkeeping.

## MultiQC Parquet Output (BETA)

Starting from version 1.29, MultiQC writes out all plot and table data in a standardized Apache Parquet file format (`multiqc.parquet`) in the `multiqc_data` directory. This feature provides several significant benefits:

- **Persistence**: The parquet file contains all the data necessary to regenerate MultiQC reports without needing access to the original analysis files
- **Reusability**: The data is structured in a way that's optimized for cross-run analysis and data warehousing
- **Interoperability**: Parquet is a widely supported columnar format that can be used with various data analysis tools and platforms

:::note
Note that the format is unstable as of 1.29 may change in 1.30, where it will be finally renamed to `multiqc.parquet`.
:::

### Parquet File Structure

The `multiqc.parquet` file contains several different types of rows that can be distinguished by the `type` column:

1. **`run_metadata`**: Contains metadata about the MultiQC run, including:
   - `creation_date`: Timestamp when the report was generated
   - `modules`: JSON-encoded list of modules included in the report
   - `data_sources`: JSON-encoded information about the data source files
   - `config`: JSON-encoded MultiQC configuration used for the run
   - `multiqc_version`: The version of MultiQC used

2. **`plot_input`**: Contains the serialized plot configuration and data:
   - `anchor`: Unique identifier for the plot
   - `plot_type`: Type of plot (e.g., "line", "bar", "heatmap", "violin", "scatter", "table")
   - `plot_input_data`: JSON-encoded representation of the plot data and configuration

3. **`table_row`**: Contains tabular data for samples and metrics:
   - `sample_name`: Name of the sample
   - `metric_name`: Name of the metric
   - `val_raw`: Raw value of the metric (numeric)
   - `val_raw_type`: Type of the raw value (e.g., "int", "float", "bool")
   - `val_str`: String representation of the value
   - `metric_col_name`: Column name in the source table
   - `module`: Name of the module that generated this data
   - `section`: Section within the module

Additional columns may be present depending on the specific plot or table type.

#### Rows and Schema

The schema is dynamically created based on the data, but here's a representative schema of the core columns:

```python
{
    "anchor": pl.Utf8,
    "type": pl.Utf8,
    "creation_date": pl.Datetime(time_unit="us"),  # no timezone specifier, but assumed UTC (for compatibility with Iceberg)
    "plot_type": pl.Utf8,
    "plot_input_data": pl.Utf8,
    "sample_name": pl.Utf8,
    "metric_name": pl.Utf8,
    "val_raw": pl.Float64,
    "val_raw_type": pl.Utf8,
    "val_str": pl.Utf8,
    "module": pl.Utf8,
    "section": pl.Utf8,
}
```

#### Working with Parquet Data

To explore the structure programmatically:

```python
import polars as pl

# Load the parquet file
df = pl.read_parquet("multiqc_data/multiqc.parquet")

# Get unique row types
print(df.select("type").unique())

# Access metadata
metadata_rows = df.filter(pl.col("type") == "run_metadata")

# Get all plot configurations
plot_inputs = df.filter(pl.col("type") == "plot_input")

# Extract tabular data from a specific module
module_data = df.filter(
    (pl.col("type") == "table_row") &
    (pl.col("module") == "fastqc")
)

# Get all unique metrics available
metrics = df.filter(pl.col("type") == "table_row").select("metric_name").unique()
```

#### Relationships Between Data

- The `anchor` column connects `plot_input` rows with their corresponding data rows
- The `module` and `section` columns in tabular data connect rows to their source modules
- `creation_date` allows tracking when the data was generated

Developers can use these relationships to reconstruct the full structure of the MultiQC report from the parquet file alone.

### Rerunning MultiQC from Parquet

One of the key benefits of the parquet output is the ability to regenerate MultiQC reports without needing the original data files:

```bash
multiqc multiqc_data/multiqc.parquet
```

This will load all the data from the parquet file and generate a new report.

Note that to be discovered, the file name must end with `*multiqc.parquet`.

### Combining Multiple MultiQC Runs

The parquet output enables easy aggregation of data from multiple MultiQC runs:

```bash
# Run MultiQC on the first set of data
multiqc /path/to/analysis1/ -o run1_output

# Run MultiQC on both the second set of data and the parquet from the first run
multiqc /path/to/analysis2/ run1_output/multiqc_data/multiqc.parquet -o combined_output
```

This will generate a report containing data from both runs. You can combine any number of parquet files with new data in a single command.

### Using MultiQC Data in Python Scripts

For programmatic access to MultiQC data, you can use the Python API to load parquet files directly:

```python
import multiqc

# Load data from a parquet file
multiqc.parse_logs('multiqc_data/multiqc.parquet')

# List loaded modules and access data
modules = multiqc.list_modules()
plots = multiqc.list_plots()
data = multiqc.get_module_data(module="fastp")
```

### Integrating with OLAP Databases

The structured format of MultiQC's parquet output makes it ideal for integration with analytical databases and OLAP systems like Apache Iceberg:

```python
import polars as pl
from pyiceberg.catalog import load_catalog

# Load the MultiQC parquet file
multiqc_df = pl.read_parquet("multiqc_data/multiqc.parquet")

# Configure and load Iceberg catalog
catalog = load_catalog(
    "glue",
    **{
        "type": "glue",
        "warehouse": "s3://your-bucket/iceberg-warehouse/"
    }
)

# Create or load Iceberg table
table = catalog.load_table("your_database.multiqc_data")

# Append data to Iceberg table
table.append(multiqc_df.to_arrow())
```

This approach enables more sophisticated analysis workflows, better reproducibility, and easier collaboration across teams - all while maintaining the comprehensive and intuitive reporting that MultiQC is known for.

### Parquet Format Options

Currently MultiQC offers two format options for the parquet output, but we might settle with only one format in the future.

1. **Long format** (default): Data is stored with columns 'sample_name', 'metric_name', 'val_raw', 'val_raw_type', and 'val_str'. This format is very flexible and ensures all data types can be preserved.

2. **Wide format**: Data is stored with each metric as a separate column, prefixed with the table name and optional namespace. While more intuitive for analytics, it may hit limits on the maximum number of columns in certain edge cases, and can have issues with mixed types (since Parquet requires columns to have consistent types).

You can configure the format in your MultiQC configuration file:

```yaml
parquet_format: "long" # or "wide"
```
