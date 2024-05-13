---
title: Using MultiQC interactive environments
description: Interactive log parsing and plotting with MultiQC
---

# Using MultiQC in interactive environments

Even though the primary way to run MultiQC is as a command line, it can also be imported
like a Python module in order to build the report interactively,
such as in a Jupyter notebook environment
(See an [example notebook](https://deploy-preview-94--multiqc.netlify.app/example-reports/jupyter/)).

MultiQC provides a set of commands to iteratively parse logs and add sections to a report.

## Functions

`multiqc.add_custom_content_section(name, anchor, description='', content_before_plot='', plot: Union[multiqc.plots.plotly.plot.Plot, str, ForwardRef(None)] = None, content='', comment='', helptext='')`
: Add a custom content section to the report. This can be used to add a custom table or other content.

    @param name: Desired section name
    @param anchor: Desired section anchor (should be unique in the session)
    @param description: Section text description
    @param content_before_plot: Content to show before the plot
    @param plot: Plot object or plot ID to show
    @param content: Content to show after the plot
    @param comment: Comment to show in the report
    @param helptext: Longer help text to show in the report, will be hidden by default, and expandable by user

`multiqc.get_general_stats_data(sample: Optional[str] = None) ‑> Dict`
: Return parsed general stats data, indexed by sample, then by data key. If sample is specified,
return only data for that sample.

    @param sample: Sample name
    @return: Dict of general stats data indexed by sample and data key

`multiqc.get_module_data(module: Optional[str] = None, sample: Optional[str] = None, key: Optional[str] = None) ‑> Dict`
: Return parsed module data, indexed (if available) by data key, then by sample. Module is either
the module name, or the anchor.

    Takes data from report.saved_raw_data, which populated by self.write_data_file() calls in modules.
    This data is not necessarily normalized, e.g. numbers can be strings or numbers, depends on
    individual module behaviour.

    @param module: Module name or anchor
    @param sample: Sample name
    @param key: Data key
    @return: Dict of module data indexed by sample and data key

`multiqc.list_data_sources() ‑> List[str]`
: Return a list of the data sources that have been loaded.

    @return: List of data sources paths from loaded modules

`multiqc.list_modules() ‑> List[str]`
: Return a list of the modules that have been loaded, in order according to config.

    @return: List of loaded module names

`multiqc.list_plots() ‑> Dict[str, List[Union[str, Dict[str, str]]]]`
: Return plot names that have been loaded, indexed by module and section.

    @return: Dict of plot names indexed by module and section

`multiqc.list_samples() ‑> List[str]`
: Return a list of the samples that have been loaded.

    @return: List of sample names from loaded modules

`multiqc.load_config(config_file: Union[str, pathlib.Path])`
: Load config on top of the current config from a MultiQC config file.

    @param config_file: Path to the config file

`multiqc.parse_data_json(path: Union[str, pathlib.Path])`
: Try find multiqc_data.json in the given directory, and load it into the report.

    @param path: Path to the directory containing multiqc_data.json or the path to the file itself.

`multiqc.parse_logs(*analysis_dir: str, verbose: Optional[bool] = None, file_list: Optional[bool] = None, prepend_dirs: Optional[bool] = None, dirs_depth: Optional[int] = None, fn_clean_sample_names: Optional[bool] = None, require_logs: Optional[bool] = None, use_filename_as_sample_name: Optional[bool] = None, strict: Optional[bool] = None, quiet: Optional[bool] = None, no_ansi: Optional[bool] = None, profile_runtime: Optional[bool] = None, no_version_check: Optional[bool] = None, ignore: List[str] = (), ignore_samples: List[str] = (), run_modules: List[str] = (), exclude_modules: List[str] = (), config_files: List[str] = (), module_order: List[Union[str, Dict]] = (), extra_fn_clean_exts: List = (), extra_fn_clean_trim: List = (), preserve_module_raw_data: bool = True)`
: Find files that MultiQC recognizes in `analysis_dir` and parse them, without generating a report.
Data can be accessed with other methods: `list_modules`, `show_plot`, `get_summarized_data`, etc.

    @param analysis_dir: Paths to search for files to parse
    @param verbose: Print more information to the console
    @param file_list: Supply a file containing a list of file paths to be searched, one per row
    @param prepend_dirs: Prepend directory to sample names
    @param dirs_depth: Prepend n directories to sample names. Negative number to take from start of path
    @param fn_clean_sample_names: Do not clean the sample names (leave as full file name)
    @param require_logs: Require all explicitly requested modules to have log files. If not, MultiQC will exit with an error
    @param use_filename_as_sample_name: Use the log filename as the sample name
    @param strict: Don't catch exceptions, run additional code checks to help development
    @param quiet: Only show log warnings
    @param no_ansi: Disable coloured log output
    @param profile_runtime: Add analysis of how long MultiQC takes to run to the report
    @param no_version_check: Disable checking the latest MultiQC version on the server
    @param ignore: Ignore analysis files
    @param ignore_samples: Ignore sample names
    @param run_modules: Use only this module. Can specify multiple times
    @param exclude_modules: Do not use this module. Can specify multiple times
    @param config_files: Specific config file to load, after those in MultiQC dir / home dir / working dir
    @param module_order: Names of modules in order of precedence to show in report
    @param extra_fn_clean_exts: Extra file extensions to clean from sample names
    @param extra_fn_clean_trim: Extra strings to clean from sample names
    @param preserve_module_raw_data: Preserve raw data from modules in the report - besides plots. Useful to use
     later interactively. Defaults to `True`. Set to `False` to save memory.

`multiqc.reset()`
: Reset the report to start fresh. Drops all previously parsed data.

`multiqc.show_plot(module: str, section: str, dataset: Optional[str] = None, flat=False, **kwargs)`
: Show a plot in the notebook.

    @param module: Module name or anchor
    @param section: Section name or anchor
    @param dataset: Dataset label, in case if plot has several tabs
    @param flat: Show plot as static images without any interactivity

`multiqc.write_report(title: Optional[str] = None, report_comment: Optional[str] = None, template: Optional[str] = None, output_dir: Optional[str] = None, filename: Optional[str] = None, make_data_dir: Optional[bool] = None, data_format: Optional[str] = None, zip_data_dir: Optional[bool] = None, force: Optional[bool] = None, make_report: Optional[bool] = None, export_plots: Optional[bool] = None, plots_force_flat: Optional[bool] = None, plots_force_interactive: Optional[bool] = None, strict: Optional[bool] = None, development: Optional[bool] = None, make_pdf: Optional[bool] = None, no_megaqc_upload: Optional[bool] = None, quiet: Optional[bool] = None, verbose: Optional[bool] = None, no_ansi: Optional[bool] = None, profile_runtime: Optional[bool] = None, no_version_check: Optional[bool] = None, run_modules: List[str] = (), exclude_modules: List[str] = (), config_files: List[str] = (), custom_css_files: List[str] = (), module_order: List[Union[str, Dict]] = ())`
: Render HTML from parsed module data, and write a report and data files to disk.

    @param title: Report title. Printed as page header, used for filename if not otherwise specified
    @param report_comment: Custom comment, will be printed at the top of the report
    @param template: Report template to use
    @param output_dir: Create report in the specified output directory
    @param filename: Report filename. Use 'stdout' to print to standard out
    @param make_data_dir: Force the parsed data directory to be created
    @param data_format: Output parsed data in a different format
    @param zip_data_dir: Compress the data directory
    @param force: Overwrite existing report and data directory
    @param make_report: Generate the report HTML. Defaults to `True`, set to `False` to only export data and plots
    @param export_plots: Export plots as static images in addition to the report
    @param plots_force_flat: Use only flat plots (static images)
    @param plots_force_interactive: Use only interactive plots (in-browser Javascript)
    @param strict: Don't catch exceptions, run additional code checks to help development
    @param development: Development mode. Do not compress and minimise JS, export uncompressed plot data
    @param make_pdf: Create PDF report. Requires Pandoc to be installed
    @param no_megaqc_upload: Don't upload generated report to MegaQC, even if MegaQC options are found
    @param quiet: Only show log warnings
    @param verbose: Print more information to the console
    @param no_ansi: Disable coloured log output
    @param profile_runtime: Add analysis of how long MultiQC takes to run to the report
    @param no_version_check: Disable checking the latest MultiQC version on the server
    @param run_modules: Use only these modules
    @param exclude_modules: Do not use these modules
    @param config_files: Specific config file to load, after those in MultiQC dir / home dir / working dir
    @param custom_css_files: Custom CSS files to include in the report
    @param module_order: Names of modules in order of precedence to show in report
