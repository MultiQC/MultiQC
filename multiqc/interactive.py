"""
This module provides functions useful to interact with MultiQC in an interactive
Python environment, such as Jupyter notebooks.
"""

import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Union, List, Optional


from multiqc import report, config
from multiqc.base_module import BaseMultiqcModule
from multiqc.core.update_config import update_config, ClConfig
from multiqc.core.file_search import file_search
from multiqc.core.exec_modules import exec_modules
from multiqc.core.version_check import check_version
from multiqc.core.write_results import write_results
from multiqc.core.exceptions import RunError
from multiqc.plots.plotly.bar import BarPlot
from multiqc.plots.plotly.box import BoxPlot
from multiqc.plots.plotly.heatmap import HeatmapPlot
from multiqc.plots.plotly.line import LinePlot
from multiqc.plots.plotly.plot import PlotType, Plot
from multiqc.plots.plotly.scatter import ScatterPlot
from multiqc.plots.plotly.violin import ViolinPlot

logger = logging.getLogger("multiqc")


def parse_logs(
    *analysis_dir: Union[str, Path],
    verbose: Optional[bool] = None,
    file_list: Optional[bool] = None,
    prepend_dirs: Optional[bool] = None,
    dirs_depth: Optional[int] = None,
    fn_clean_sample_names: Optional[bool] = None,
    require_logs: Optional[bool] = None,
    use_filename_as_sample_name: Optional[bool] = None,
    strict: Optional[bool] = None,
    quiet: Optional[bool] = None,
    no_ansi: Optional[bool] = None,
    profile_runtime: Optional[bool] = None,
    no_version_check: Optional[bool] = None,
    ignore: List[str] = (),
    ignore_samples: List[str] = (),
    run_modules: List[str] = (),
    exclude_modules: List[str] = (),
    config_files: List[str] = (),
    module_order: List[Union[str, Dict]] = (),
    extra_fn_clean_exts: List = (),
    extra_fn_clean_trim: List = (),
    preserve_module_raw_data: bool = True,
):
    """
    Find files that MultiQC recognizes in `analysis_dir` and parse them, without generating a report.
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
    """
    assert isinstance(analysis_dir, tuple)
    if not all(isinstance(d, (str, Path)) for d in analysis_dir):
        raise ValueError("Path arguments should be path-like or strings, got:", analysis_dir)

    update_config(*analysis_dir, cfg=ClConfig(**{k: v for k, v in locals().items() if k != "analysis_dir"}))

    check_version(parse_logs.__name__)

    report.reset_file_search()
    try:
        searched_modules = file_search()
        exec_modules(searched_modules, clean_up=False)
    except RunError as e:
        if e.message:
            logger.critical(e.message)


def parse_data_json(path: Union[str, Path]):
    """
    Try find multiqc_data.json in the given directory, and load it into the report.

    @param path: Path to the directory containing multiqc_data.json or the path to the file itself.
    """
    check_version(parse_data_json.__name__)

    json_path_found = False
    json_path: Path
    if path.endswith(".json"):
        json_path = path
        json_path_found = True
    else:
        json_path = Path(path) / "multiqc_data.json"
        if json_path.exists():
            json_path_found = True

    if not json_path_found:
        logger.error(f"multiqc_data.json not found in {path}")
        return

    logger.info(f"Loading data from {json_path}")
    try:
        with json_path.open("r") as f:
            data = json.load(f)

        for mod, sections in data["report_data_sources"].items():
            logger.info(f"Loaded module {mod}")
            for section, sources in sections.items():
                for sname, source in sources.items():
                    report.data_sources[mod][section][sname] = source
        for id, plot_dump in data["report_plot_data"].items():
            logger.info(f"Loaded plot {id}")
            report.plot_data[id] = plot_dump
    except (json.JSONDecodeError, KeyError) as e:
        logger.error(f"Error loading data from multiqc_data.json: {e}")


def list_data_sources() -> List[str]:
    """
    Return a list of the data sources that have been loaded.

    @return: List of data sources paths from loaded modules
    """
    file_list = []
    for mod, sections in report.data_sources.items():
        for section, sources in sections.items():
            for sname, source in sources.items():
                file_list.append(source)
    return file_list


def list_modules() -> List[str]:
    """
    Return a list of the modules that have been loaded, in order according to config.

    @return: List of loaded module names
    """
    return [m.name for m in report.modules]


def list_samples() -> List[str]:
    """
    Return a list of the samples that have been loaded.

    @return: List of sample names from loaded modules
    """
    samples = set()

    for mod, sections in report.data_sources.items():
        for section, sources in sections.items():
            for sname, source in sources.items():
                samples.add(sname)

    return sorted(samples)


def list_plots() -> Dict[str, List[Union[str, Dict[str, str]]]]:
    """
    Return plot names that have been loaded, indexed by module and section.

    @return: Dict of plot names indexed by module and section
    """

    result = dict()
    for module in report.modules:
        result[module.name]: List[Union[str, Dict[str, str]]] = list()
        for section in module.sections:
            if not section.plot_id:
                continue
            section_id = section.name or section.anchor
            plot_id = section.plot_id
            if plot_id not in report.plot_by_id:
                raise ValueError(f'CRITICAL: Plot "{plot_id}" not found in report.plot_by_id')
            plot = report.plot_by_id[plot_id]
            if len(plot.datasets) == 1:
                result[module.name].append(section_id)
            if len(plot.datasets) > 1:
                result[module.name].append({section_id: [d.label for d in plot.datasets]})

    return result


def get_plot(
    module: str,
    section: str,
) -> Plot:
    """
    Get plot Object by module name and section ID.

    @param module: Module name or anchor
    @param section: Section name or anchor
    """
    mod = next((m for m in report.modules if m.name == module or m.anchor == module), None)
    if not mod:
        raise ValueError(f'Module "{module}" is not found. Use multiqc.list_modules() to list available modules')

    sec = next((s for s in mod.sections if (s.name and s.name == section) or s.anchor == section), None)
    if not sec:
        raise ValueError(f'Section "{section}" is not found in module "{module}"')

    if sec.plot_id is None:
        raise ValueError(f"Section {section} doesn't contain a Plot object")

    return report.plot_by_id[sec.plot_id]


def _load_plot(dump: Dict) -> Plot:
    """
    Load a plot and datasets from a JSON dump.
    """

    plot_type = PlotType(dump["plot_type"])
    if plot_type == PlotType.LINE:
        return LinePlot(**dump)
    elif plot_type == PlotType.BAR:
        return BarPlot(**dump)
    elif plot_type == PlotType.BOX:
        return BoxPlot(**dump)
    elif plot_type == PlotType.SCATTER:
        return ScatterPlot(**dump)
    elif plot_type == PlotType.HEATMAP:
        return HeatmapPlot(**dump)
    elif plot_type == PlotType.VIOLIN:
        return ViolinPlot(**dump)
    else:
        raise ValueError(f"Plot type {plot_type} is unknown or unsupported")


def get_general_stats_data(sample: Optional[str] = None) -> Dict:
    """
    Return parsed general stats data, indexed by sample, then by data key. If sample is specified,
    return only data for that sample.

    @param sample: Sample name
    @return: Dict of general stats data indexed by sample and data key
    """

    data = defaultdict(dict)
    for data_by_sample, header in zip(report.general_stats_data, report.general_stats_headers):
        for s, val_by_key in data_by_sample.items():
            if sample and s != sample:
                continue
            for key, val in val_by_key.items():
                if key in header:
                    key = f"{header[key].get('namespace', '')}.{key}"
                    data[s][key] = val
    if sample:
        if not data:
            return {}
        return data[sample]

    return data


def get_module_data(
    module: Optional[str] = None,
    sample: Optional[str] = None,
    key: Optional[str] = None,
) -> Dict:
    """
    Return parsed module data, indexed (if available) by data key, then by sample. Module is either
    the module name, or the anchor.

    Takes data from report.saved_raw_data, which populated by self.write_data_file() calls in modules.
    This data is not necessarily normalized, e.g. numbers can be strings or numbers, depends on
    individual module behaviour.

    @param module: Module name or anchor
    @param sample: Sample name
    @param key: Data key
    @return: Dict of module data indexed by sample and data key
    """

    if sample and sample not in list_samples():
        raise ValueError(f"Sample '{sample}' is not found. Use multiqc.list_samples() to list available samples")
    if module and module not in list_modules():
        raise ValueError(f"Module '{module}' is not found. Use multiqc.list_modules() to list available modules")

    data_by_module = {}
    for m in report.modules:
        if module and (m.name != module and m.anchor != module):
            continue

        module_data = m.saved_raw_data
        if sample:
            module_data = {k: v.get(sample, {}) for k, v in module_data.items()}
        if key:
            if module and key not in m.saved_raw_data:
                raise ValueError(f"Key '{key}' is not found in module '{module}'")
            module_data = module_data.get(key, {})
        elif len(module_data) == 1:  # only one key, flatten
            module_data = module_data[list(module_data.keys())[0]]

        data_by_module[m.name] = module_data

    if module:
        if not data_by_module:
            return {}
        return data_by_module[module]

    return data_by_module


def reset():
    """
    Reset the report to start fresh. Drops all previously parsed data.
    """

    config.reset()
    report.reset()


def add_custom_content_section(
    name,
    anchor,
    description="",
    content_before_plot="",
    plot: Optional[Union[Plot, str]] = None,
    content="",
    comment="",
    helptext="",
):
    """
    Add a custom content section to the report. This can be used to add a custom table or other content.

    @param name: Desired section name
    @param anchor: Desired section anchor (should be unique in the session)
    @param description: Section text description
    @param content_before_plot: Content to show before the plot
    @param plot: Plot object or plot ID to show
    @param content: Content to show after the plot
    @param comment: Comment to show in the report
    @param helptext: Longer help text to show in the report, will be hidden by default, and expandable by user
    """

    module = BaseMultiqcModule(
        name=name,
        anchor=anchor,
        info=description,
        comment=comment,
    )
    module.add_section(
        name=name,
        anchor=anchor,
        description=description,
        helptext=helptext,
        content_before_plot=content_before_plot,
        plot=plot,
        content=content,
        comment=comment,
    )
    report.modules.append(module)


def write_report(
    title: Optional[str] = None,
    report_comment: Optional[str] = None,
    template: Optional[str] = None,
    output_dir: Optional[str] = None,
    filename: Optional[str] = None,
    make_data_dir: Optional[bool] = None,
    data_format: Optional[str] = None,
    zip_data_dir: Optional[bool] = None,
    force: Optional[bool] = None,
    overwrite: Optional[bool] = None,
    make_report: Optional[bool] = None,
    export_plots: Optional[bool] = None,
    plots_force_flat: Optional[bool] = None,
    plots_force_interactive: Optional[bool] = None,
    strict: Optional[bool] = None,
    development: Optional[bool] = None,
    make_pdf: Optional[bool] = None,
    no_megaqc_upload: Optional[bool] = None,
    quiet: Optional[bool] = None,
    verbose: Optional[bool] = None,
    no_ansi: Optional[bool] = None,
    profile_runtime: Optional[bool] = None,
    no_version_check: Optional[bool] = None,
    run_modules: List[str] = (),
    exclude_modules: List[str] = (),
    config_files: List[str] = (),
    custom_css_files: List[str] = (),
    module_order: List[Union[str, Dict]] = (),
):
    """
    Render HTML from parsed module data, and write a report and data files to disk.

    @param title: Report title. Printed as page header, used for filename if not otherwise specified
    @param report_comment: Custom comment, will be printed at the top of the report
    @param template: Report template to use
    @param output_dir: Create report in the specified output directory
    @param filename: Report filename. Use 'stdout' to print to standard out
    @param make_data_dir: Force the parsed data directory to be created
    @param data_format: Output parsed data in a different format
    @param zip_data_dir: Compress the data directory
    @param force: Overwrite existing report and data directory
    @param overwrite: Same as force
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
    """

    if force is None and overwrite is not None:
        force = overwrite
    params = locals()
    del params["overwrite"]
    update_config(cfg=ClConfig(**params))

    check_version(write_report.__name__)

    if len(report.modules) == 0:
        logger.error("No analysis results found to make a report")
        return

    try:
        write_results()

    except RunError as e:
        if e.message:
            logger.critical(e.message)


def load_config(config_file: Union[str, Path]):
    """
    Load config on top of the current config from a MultiQC config file.

    @param config_file: Path to the config file
    """

    update_config()

    path = Path(config_file)
    if not path.exists():
        raise ValueError(f"Config file '{config_file}' not found")

    config.session_user_config_files.append(path.absolute())
    config.load_config_file(config_file)
