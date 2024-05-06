---
title: Using MultiQC interactive environments
description: Interactive log parsing and plotting with MultiQC
---

# Using MultiQC in interactive environments

Even though the primary way to run MultiQC is as a command line, it can also be imported
like a Python module in order to build the report interactively,
e.g. in a [Jupyter notebook environment](https://deploy-preview-94--multiqc.netlify.app/example-reports/jupyter/).

MultiQC provides a set of commands to iteratively parse logs and add sections to a report.

`add_custom_content_section(name, anchor, description='', content_before_plot='', plot: Union[multiqc.plots.plotly.plot.Plot, str, ForwardRef(None)] = None, content='', comment='', helptext='')`
: Add a custom content section to the report. This can be used to add a custom table or other content.

`get_general_stats_data(sample: Optional[str] = None) ‑> Dict`
: Return parsed general stats data, indexed by sample, then by data key. If sample is specified,
return only data for that sample.

`get_module_data(module: Optional[str] = None, sample: Optional[str] = None, key: Optional[str] = None) ‑> Dict`
: Return parsed module data, indexed (if available) by data key, then by sample. Module is either
the module name, or the anchor.

    Takes data from report.saved_raw_data, which populated by self.write_data_file() calls in modules.
    This data is not necessarily normalized, e.g. numbers can be strings or numbers, depends on
    individual module behaviour.

`list_data_sources() ‑> List[str]`
: Return a list of the data sources that have been loaded.

`list_modules() ‑> List[str]`
: Return a list of the modules that have been loaded, in order according to config.

`list_plots() ‑> Dict[str, List[Union[str, Dict[str, str]]]]`
: Return plot names that have been loaded, indexed by module and section.

`list_samples() ‑> List[str]`
: Return a list of the samples that have been loaded.

`load_config(config_file: Union[str, pathlib.Path])`
: Load config on top of the current config from a MultiQC config file.

`parse_data_json(path: Union[str, pathlib.Path])`
: Try find multiqc_data.json in the given directory, and load it into the report.

`parse_logs(*analysis_dir: str, verbose: Optional[bool] = None, file_list: Optional[bool] = None, prepend_dirs: Optional[bool] = None, dirs_depth: Optional[int] = None, fn_clean_sample_names: Optional[bool] = None, require_logs: Optional[bool] = None, use_filename_as_sample_name: Optional[bool] = None, strict: Optional[bool] = None, quiet: Optional[bool] = None, no_ansi: Optional[bool] = None, profile_runtime: Optional[bool] = None, no_version_check: Optional[bool] = None, ignore: List[str] = (), ignore_samples: List[str] = (), run_modules: List[str] = (), exclude_modules: List[str] = (), config_files: List[str] = (), module_order: List[Union[str, Dict]] = (), extra_fn_clean_exts: List = (), extra_fn_clean_trim: List = ())`
: Find files that MultiQC recognizes in `analysis_dir` and parse them, without generating a report.
Data can be accessed with other methods: `list_modules`, `show_plot`, `get_summarized_data`, etc.

`reset()`
: Reset the report to start fresh. Drops all previously parsed data.

`show_plot(module: str, section: str, dataset: Optional[str] = None, flat=False, **kwargs)`
: Show a plot in the notebook.

`write_report(title: Optional[str] = None, report_comment: Optional[str] = None, template: Optional[str] = None, output_dir: Optional[str] = None, filename: Optional[str] = None, make_data_dir: Optional[bool] = None, data_format: Optional[str] = None, zip_data_dir: Optional[bool] = None, force: Optional[bool] = None, make_report: Optional[bool] = None, export_plots: Optional[bool] = None, plots_force_flat: Optional[bool] = None, plots_force_interactive: Optional[bool] = None, strict: Optional[bool] = None, development: Optional[bool] = None, make_pdf: Optional[bool] = None, no_megaqc_upload: Optional[bool] = None, quiet: Optional[bool] = None, verbose: Optional[bool] = None, no_ansi: Optional[bool] = None, profile_runtime: Optional[bool] = None, no_version_check: Optional[bool] = None, run_modules: List[str] = (), exclude_modules: List[str] = (), config_files: List[str] = (), custom_css_files: List[str] = (), module_order: List[Union[str, Dict]] = ())`
: Render HTML from parsed module data, and write a report and data files to disk.
