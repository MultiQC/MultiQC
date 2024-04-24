"""
multiqc.multiqc
~~~~~~~~~~~~~~~~~~~~~
The main function to run MultiQC. Sorry about the messy namespace.
Primarily called by multiqc.__main__.py
Imported by __init__.py so available as multiqc.run()
"""

import json
import logging
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Union, List, Optional

import rich_click as click

from multiqc.core.write_results import _write_html_and_data, _write_results
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots.plotly.bar import BarPlot
from multiqc.plots.plotly.box import BoxPlot
from multiqc.plots.plotly.heatmap import HeatmapPlot
from multiqc.plots.plotly.line import LinePlot
from multiqc.plots.plotly.plot import go, PlotType, Plot
from multiqc.plots.plotly.scatter import ScatterPlot
from multiqc.plots.plotly.violin import ViolinPlot
from multiqc.utils import config, plugin_hooks, report, util_functions, log
from multiqc.core.cl_to_config import _cl_to_config
from multiqc.core.file_search import _file_search
from multiqc.core.run_modules import _run_modules
from multiqc.core.utils import RunResult, _RunError

# Set up logging
start_execution_time = time.time()
logger = logging.getLogger(__name__)


# Configuration for rich-click CLI help
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.HEADER_TEXT = f"[dark_orange]///[/] [bold][link=https://multiqc.info]MultiQC[/link][/] :{util_functions.choose_emoji()}: [dim]| v{config.version}"
click.rich_click.FOOTER_TEXT = "See [link=http://multiqc.info]http://multiqc.info[/] for more details."
click.rich_click.ERRORS_SUGGESTION = f"This is MultiQC [cyan]v{config.version}[/]\nFor more help, run '[yellow]multiqc --help[/]' or visit [link=http://multiqc.info]http://multiqc.info[/]"
click.rich_click.STYLE_ERRORS_SUGGESTION = ""
click.rich_click.OPTION_GROUPS = {
    "multiqc": [
        {
            "name": "Main options",
            "options": [
                "--force",
                "--config",
                "--cl-config",
                "--filename",
                "--outdir",
                "--ignore",
                "--ignore-samples",
                "--ignore-symlinks",
                "--file-list",
            ],
        },
        {
            "name": "Choosing modules to run",
            "options": [
                "--module",
                "--exclude",
            ],
        },
        {
            "name": "Sample handling",
            "options": [
                "--dirs",
                "--dirs-depth",
                "--fullnames",
                "--fn_as_s_name",
                "--rename-sample-names",
                "--replace-names",
            ],
        },
        {
            "name": "Report customisation",
            "options": [
                "--title",
                "--comment",
                "--template",
                "--sample-names",
                "--sample-filters",
                "--custom-css-file",
            ],
        },
        {
            "name": "Output files",
            "options": [
                "--flat",
                "--interactive",
                "--export",
                "--data-dir",
                "--no-data-dir",
                "--data-format",
                "--zip-data-dir",
                "--no-report",
                "--pdf",
            ],
        },
        {
            "name": "MultiQC behaviour",
            "options": [
                "--verbose",
                "--quiet",
                "--strict",
                "--development",
                "--require-logs",
                "--profile-runtime",
                "--no-megaqc-upload",
                "--no-ansi",
                "--version",
                "--help",
            ],
        },
    ],
}


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument("analysis_dir", type=click.Path(exists=True), nargs=-1, required=True, metavar="[ANALYSIS DIRECTORY]")
@click.option("-f", "--force", is_flag=True, help="Overwrite any existing reports")
@click.option("-d", "--dirs", "prepend_dirs", is_flag=True, help="Prepend directory to sample names")
@click.option(
    "-dd",
    "--dirs-depth",
    "dirs_depth",
    type=int,
    help="Prepend [yellow i]n[/] directories to sample names. Negative number to take from start of path.",
)
@click.option("-s", "" "transformation", flag_value="upper", default=True)
@click.option(
    "-s",
    "--fullnames",
    "fn_clean_sample_names",
    flag_value=False,
    help="Do not clean the sample names [i](leave as full file name)[/]",
)
@click.option(
    "-i",
    "--title",
    type=str,
    help="Report title. Printed as page header, used for filename if not otherwise specified.",
)
@click.option(
    "-b", "--comment", "report_comment", type=str, help="Custom comment, will be printed at the top of the report."
)
@click.option("-n", "--filename", type=str, help="Report filename. Use '[yellow i]stdout[/]' to print to standard out.")
@click.option("-o", "--outdir", type=str, help="Create report in the specified output directory.")
@click.option(
    "-t", "--template", type=click.Choice(config.avail_templates), metavar=None, help="Report template to use."
)
@click.option("-x", "--ignore", type=str, multiple=True, metavar="GLOB EXPRESSION", help="Ignore analysis files")
@click.option(
    "--ignore-samples", "ignore_samples", type=str, multiple=True, metavar="GLOB EXPRESSION", help="Ignore sample names"
)
@click.option("--ignore-symlinks", "ignore_symlinks", is_flag=True, help="Ignore symlinked directories and files")
@click.option(
    "--fn_as_s_name", "use_filename_as_sample_name", is_flag=True, help="Use the log filename as the sample name"
)
@click.option(
    "--replace-names",
    "replace_names",
    type=click.Path(exists=True, readable=True),
    help="TSV file to rename sample names during report generation",
)
@click.option(
    "--sample-names",
    "sample_names",
    type=click.Path(exists=True, readable=True),
    help="TSV file containing alternative sample names for renaming buttons in the report",
)
@click.option(
    "--sample-filters",
    "sample_filters",
    type=click.Path(exists=True, readable=True),
    help="TSV file containing show/hide patterns for the report",
)
@click.option(
    "-l", "--file-list", is_flag=True, help="Supply a file containing a list of file paths to be searched, one per row"
)
@click.option(
    "-e",
    "--exclude",
    metavar="[MODULE NAME]",
    type=click.Choice(sorted(["general_stats"] + list(config.avail_modules.keys()))),
    multiple=True,
    help="Do not use this module. Can specify multiple times.",
)
@click.option(
    "-m",
    "--module",
    metavar="[MODULE NAME]",
    type=click.Choice(sorted(config.avail_modules.keys())),
    multiple=True,
    help="Use only this module. Can specify multiple times.",
)
@click.option(
    "--require-logs",
    "require_logs",
    is_flag=True,
    help="Require all explicitly requested modules to have log files. If not, MultiQC will exit with an error.",
)
@click.option(
    "--data-dir/--no-data-dir", "make_data_dir", is_flag=True, help="Force the parsed data directory to be created."
)
@click.option(
    "-k",
    "--data-format",
    "data_format",
    type=click.Choice(list(config.data_format_extensions.keys())),
    help="Output parsed data in a different format.",
)
@click.option("-z", "--zip-data-dir", "zip_data_dir", is_flag=True, help="Compress the data directory.")
@click.option(
    "--no-report", "make_report", flag_value=False, help="Do not generate a report, only export data and plots"
)
@click.option(
    "-p", "--export", "export_plots", is_flag=True, help="Export plots as static images in addition to the report"
)
@click.option("-fp", "--flat", "plots_flat", is_flag=True, help="Use only flat plots [i](static images)[/]")
@click.option(
    "-ip",
    "--interactive",
    "plots_interactive",
    is_flag=True,
    help="Use only interactive plots [i](in-browser Javascript)[/]",
)
@click.option(
    "--strict",
    "strict",
    is_flag=True,
    help="Don't catch exceptions, run additional code checks to help development.",
)
@click.option("--lint", "lint", is_flag=True, hidden=True, help="DEPRECATED: use --strict instead")
@click.option(
    "--development",
    "--dev",
    "development",
    is_flag=True,
    help="Development mode. Do not compress and minimise JS, export uncompressed plot data",
)
@click.option(
    "--pdf",
    "make_pdf",
    is_flag=True,
    help="Creates PDF report with the [i]'simple'[/] template. Requires [link=https://pandoc.org/]Pandoc[/] to be installed.",
)
@click.option(
    "--no-megaqc-upload",
    "no_megaqc_upload",
    is_flag=True,
    help="Don't upload generated report to MegaQC, even if MegaQC options are found",
)
@click.option(
    "-c",
    "--config",
    "config_file",
    type=click.Path(exists=True, readable=True),
    multiple=True,
    help="Specific config file to load, after those in MultiQC dir / home dir / working dir.",
)
@click.option("--cl-config", type=str, multiple=True, help="Specify MultiQC config YAML on the command line")
@click.option("-v", "--verbose", count=True, default=0, help="Increase output verbosity.")
@click.option("-q", "--quiet", is_flag=True, help="Only show log warnings")
@click.option("--profile-runtime", is_flag=True, help="Add analysis of how long MultiQC takes to run to the report")
@click.option("--no-ansi", is_flag=True, help="Disable coloured log output")
@click.option(
    "--custom-css-file",
    "custom_css_files",
    type=click.Path(exists=True, readable=True),
    multiple=True,
    help="Custom CSS file to add to the final report",
)
@click.version_option(config.version, prog_name="multiqc")
def run_cli(**kwargs):
    # Main MultiQC run command for use with the click command line, complete with all click function decorators.
    # To make it easy to use MultiQC within notebooks and other locations that don't need click, we simply pass the
    # parsed variables on to a vanilla python function.

    """MultiQC aggregates results from bioinformatics analyses across many samples into a single report.

    It searches a given directory for analysis logs and compiles a HTML report.
    It's a general use tool, perfect for summarising the output from numerous
    bioinformatics tools.

    To run, supply with one or more directory to scan for analysis results.
    For example, to run in the current working directory, use '[blue bold]multiqc .[/]'
    """

    # Pass on to a regular function that can be used easily without click
    multiqc_run = run(**kwargs)

    # End execution using the exit code returned from MultiQC
    sys.exit(multiqc_run.sys_exit_code)


def parse_logs(analysis_dir, *args, **kwargs) -> RunResult:
    """
    Parse files without generating a report. Useful to work with MultiQC interactively. Data can be accessed
    with other methods: `list_modules`, `show_plot`, `get_summarized_data`, etc.
    """
    # First, try find multiqc_data.json in the given directory and load it into the report.
    res = _load_multiqc_data_json(analysis_dir, *args, **kwargs)
    if res:
        return res

    # multiqc_data.json was not found, so proceed to the standard run of finding logs
    # and running modules on them.
    kwargs["no_report"] = True  # Disable report generation

    return run(analysis_dir, *args, **kwargs)


def _load_multiqc_data_json(analysis_dir, *args, **kwargs) -> Optional[RunResult]:
    json_path_found = False
    json_path = None
    if not isinstance(analysis_dir, list):
        json_path = Path("multiqc_data.json")
        if json_path.exists():
            json_path_found = True
        else:
            json_path = Path(analysis_dir) / "multiqc_data.json"
            if json_path.exists():
                json_path_found = True

    if not json_path_found:
        return None

    # Loading from previous JSON
    _cl_to_config(analysis_dir=None, *args, **kwargs)
    logger.info(f"Loading data from {json_path}")
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
    return RunResult()


def list_modules() -> List[str]:
    """
    Return a list of the modules that have been loaded.
    """
    return [m.name for m in report.modules_output]


def list_data_sources() -> List[str]:
    """
    Return a list of the data sources that have been loaded.
    """
    file_list = []
    for mod, sections in report.data_sources.items():
        for section, sources in sections.items():
            for sname, source in sources.items():
                file_list.append(source)
    return file_list


def list_samples() -> List[str]:
    """
    Return a list of the samples that have been loaded.
    """
    sample_names = set()
    for mod, sections in report.data_sources.items():
        for section, sources in sections.items():
            for sname, source in sources.items():
                sample_names.add(sname)
    return sorted(sample_names)


def list_plots() -> List[str]:
    """
    Return a list of the plots that have been loaded for a given module,
    along with the number of datasets in each plot.
    """
    return list(
        f"{id} ({len(plot['datasets'])} datasets)" if len(plot["datasets"]) > 1 else id
        for id, plot in report.plot_data.items()
    )


def _load_plot(dump: Dict) -> Plot:
    """
    Load a plot and datasets from a JSON dump.
    """
    plot_type = PlotType(dump["plot_type"])
    dump["layout"] = go.Layout(dump["layout"])
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


def show_plot(plot_id: str, dataset_id=0, **kwargs) -> go.Figure:
    """
    Show a plot in the notebook.
    """
    dump = report.plot_data[plot_id]
    model = _load_plot(dump)
    return model.show(dataset_id=dataset_id, **kwargs)


def get_general_stats_data(sample: Optional[str] = None) -> Dict:
    """
    Return parsed general stats data indexed by sample, then by data key. If sample is specified, return only data
    for that sample.
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
            raise ValueError(f"Sample '{sample}' not found in loaded data. Available samples: {list_samples()}")
        return data[sample]

    return data


def get_module_data(module: Optional[str] = None, sample: Optional[str] = None, key: Optional[str] = None) -> Dict:
    """
    Return parsed module data, indexed (optionally) by data key, then by sample. Module is either the module
    name, or the anchor.
    """
    data_by_module = {}
    for m in report.modules_output:
        if module and (m.name != module and m.anchor != module):
            continue

        module_data = m.saved_raw_data
        if sample:
            module_data = {k: v.get(sample, {}) for k, v in module_data.items()}
        if key:
            module_data = module_data.get(key, {})
        elif len(module_data) == 1:  # only one key, flatten
            module_data = module_data[list(module_data.keys())[0]]

        data_by_module[m.name] = module_data

    if module:
        if not data_by_module:
            raise ValueError(
                f"Module '{module}' not found in loaded data. Available modules: {list(data_by_module.keys())}"
            )
        return data_by_module[module]

    return data_by_module


def reset():
    """
    Reset the report to start fresh. Drops all previously parsed data.
    """
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


def write_report():
    """
    Write HTML and data files to disk. Useful to work with MultiQC interactively, after loading data with `load`.
    """
    config.make_report = True
    _write_html_and_data(
        tmp_dir=report.tmp_dir,
        filename=config.filename,
    )


# Main function that runs MultiQC. Available to use within an interactive Python environment
def run(
    analysis_dir,
    prepend_dirs=None,
    dirs_depth=None,
    fn_clean_sample_names=None,
    title=None,
    report_comment=None,
    template=None,
    module=(),
    require_logs=None,
    exclude=(),
    outdir=None,
    ignore=(),
    ignore_samples=(),
    use_filename_as_sample_name=None,
    replace_names=None,
    sample_names=None,
    sample_filters=None,
    file_list=False,
    filename=None,
    make_data_dir=None,
    data_format=None,
    zip_data_dir=None,
    force=None,
    ignore_symlinks=None,
    make_report=None,
    export_plots=None,
    plots_force_flat=None,
    plots_force_interactive=None,
    strict=None,
    lint=None,  # Deprecated since v1.17
    development=None,
    make_pdf=None,
    no_megaqc_upload=None,
    config_file=(),
    cl_config=(),
    verbose=None,
    quiet=None,
    no_ansi=None,
    profile_runtime=None,
    custom_css_files=(),
    clean_up=True,
    **kwargs,
) -> RunResult:
    """
    MultiQC aggregates results from bioinformatics analyses across many samples into a single report.

    It searches a given directory for analysis logs and compiles a HTML report.
    It's a general use tool, perfect for summarising the output from numerous
    bioinformatics tools.

    To run, supply with one or more directory to scan for analysis results.
    To run here, use 'multiqc .'

    See http://multiqc.info for more details.

    Author: Phil Ewels (http://phil.ewels.co.uk)
    """

    _cl_to_config(
        analysis_dir=analysis_dir,
        prepend_dirs=prepend_dirs,
        dirs_depth=dirs_depth,
        fn_clean_sample_names=fn_clean_sample_names,
        title=title,
        report_comment=report_comment,
        template=template,
        module=module,
        require_logs=require_logs,
        exclude=exclude,
        outdir=outdir,
        use_filename_as_sample_name=use_filename_as_sample_name,
        replace_names=replace_names,
        sample_names=sample_names,
        sample_filters=sample_filters,
        filename=filename,
        make_data_dir=make_data_dir,
        data_format=data_format,
        zip_data_dir=zip_data_dir,
        force=force,
        ignore_symlinks=ignore_symlinks,
        make_report=make_report,
        export_plots=export_plots,
        plots_force_flat=plots_force_flat,
        plots_force_interactive=plots_force_interactive,
        strict=strict,
        lint=lint,
        development=development,
        make_pdf=make_pdf,
        no_megaqc_upload=no_megaqc_upload,
        config_file=config_file,
        cl_config=cl_config,
        profile_runtime=profile_runtime,
        quiet=quiet,
        verbose=verbose,
        no_ansi=no_ansi,
        custom_css_files=custom_css_files,
        custom_config=kwargs,
    )

    try:
        run_modules, run_module_names = _file_search(
            file_list=file_list,
            ignore=ignore,
            ignore_samples=ignore_samples,
        )

        _run_modules(
            tmp_dir=report.tmp_dir,
            filename=filename,
            run_modules=run_modules,
            run_module_names=run_module_names,
        )

        _write_results(
            filename=filename,
        )

    except _RunError as e:
        if e.message:
            logger.critical(e.message)
        return RunResult(message=e.message, sys_exit_code=e.sys_exit_code)

    plugin_hooks.mqc_trigger("execution_finish")

    report.runtimes["total"] = time.time() - start_execution_time
    if config.profile_runtime:
        logger.warning(f"Run took {report.runtimes['total']:.2f} seconds")
        logger.warning(f" - {report.runtimes['total_sp']:.2f}s: Searching files")
        logger.warning(f" - {report.runtimes['total_mods']:.2f}s: Running modules")
        if config.make_report:
            logger.warning(f" - {report.runtimes['total_compression']:.2f}s: Compressing report data")
            logger.info(f"For more information, see the 'Run Time' section in {os.path.relpath(config.output_fn)}")

    if report.num_flat_plots > 0 and not config.plots_force_flat:
        if not config.plots_force_interactive:
            log.rich_console.print(
                "[blue]|           multiqc[/] | "
                "Flat-image plots used. Disable with '--interactive'. "
                "See [link=https://multiqc.info/docs/#flat--interactive-plots]docs[/link]."
            )

    sys_exit_code = 0
    if config.strict and len(report.lint_errors) > 0:
        logger.error(f"Found {len(report.lint_errors)} linting errors!\n" + "\n".join(report.lint_errors))
        sys_exit_code = 1

    logger.info("MultiQC complete")

    if clean_up:
        # Move the log file into the data directory
        log.move_tmp_log()

    return RunResult(sys_exit_code=sys_exit_code)
