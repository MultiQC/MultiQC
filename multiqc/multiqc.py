"""
multiqc.multiqc
~~~~~~~~~~~~~~~~~~~~~
The main function to run MultiQC. Sorry about the messy namespace.
Primarily called by multiqc.__main__.py
Imported by __init__.py so available as multiqc.run()
"""

import base64
import errno
import io
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
import json
import traceback
from importlib.metadata import EntryPoint
from typing import Dict, Union, Callable, List, Optional, Tuple

import jinja2
import requests
import rich
import rich_click as click
from packaging import version
from rich.syntax import Syntax

from .modules.base_module import ModuleNoSamplesFound, BaseMultiqcModule
from .plots import table
from .utils import config, log, megaqc, plugin_hooks, report, software_versions, strict_helpers, util_functions
from .utils.util_functions import strtobool

# Set up logging
start_execution_time = time.time()
logger = config.logger

OLDEST_SUPPORTED_PYTHON_VERSION = "3.8"


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
@click.option("-d", "--dirs", is_flag=True, help="Prepend directory to sample names")
@click.option(
    "-dd",
    "--dirs-depth",
    "dirs_depth",
    type=int,
    help="Prepend [yellow i]n[/] directories to sample names. Negative number to take from start of path.",
)
@click.option(
    "-s",
    "--fullnames",
    "no_clean_sname",
    is_flag=True,
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
@click.option("--data-dir", "make_data_dir", is_flag=True, help="Force the parsed data directory to be created.")
@click.option(
    "--no-data-dir", "no_data_dir", is_flag=True, help="Prevent the parsed data directory from being created."
)
@click.option(
    "-k",
    "--data-format",
    "data_format",
    type=click.Choice(list(config.data_format_extensions.keys())),
    help="Output parsed data in a different format.",
)
@click.option("-z", "--zip-data-dir", "zip_data_dir", is_flag=True, help="Compress the data directory.")
@click.option("--no-report", "no_report", is_flag=True, help="Do not generate a report, only export data and plots")
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


class RunResult:
    """
    Returned by a MultiQC run for interactive use. Contains the following information:

    * appropriate error code (e.g. 1 if a module broke, 0 on success)
    * error message if a module broke
    * report instance
    * config instance

    """

    def __init__(
        self,
        sys_exit_code: int = 0,
        message: str = "",
    ):
        self.sys_exit_code = sys_exit_code
        self.message = message
        self.report = report
        self.config = config


class _RunError(Exception):
    """
    Used internally in `run` to pass errors from sub-steps.
    """

    def __init__(self, message: str = "", sys_exit_code: int = 1):
        self.message = message
        self.sys_exit_code = sys_exit_code


# Main function that runs MultiQC. Available to use within an interactive Python environment
def run(
    analysis_dir,
    dirs=False,
    dirs_depth=None,
    no_clean_sname=False,
    title=None,
    report_comment=None,
    template=None,
    module=(),
    require_logs=False,
    exclude=(),
    outdir=None,
    ignore=(),
    ignore_samples=(),
    use_filename_as_sample_name=False,
    replace_names=None,
    sample_names=None,
    sample_filters=None,
    file_list=False,
    filename=None,
    make_data_dir=False,
    no_data_dir=False,
    data_format=None,
    zip_data_dir=False,
    force=True,
    ignore_symlinks=False,
    no_report=False,
    export_plots=False,
    plots_flat=False,
    plots_interactive=False,
    strict=False,
    lint=False,  # Deprecated since v1.17
    development=False,
    make_pdf=False,
    no_megaqc_upload=False,
    config_file=(),
    cl_config=(),
    verbose=0,
    quiet=False,
    profile_runtime=False,
    no_ansi=False,
    custom_css_files=(),
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

    # Throw an error if we are using an unsupported version of Python
    if sys.version_info < tuple(map(int, OLDEST_SUPPORTED_PYTHON_VERSION.split("."))):
        return RunResult(
            message="You are running MultiQC with Python {}. "
            "Please upgrade Python! MultiQC does not support Python < {}, "
            "things will break.".format(sys.version_info, OLDEST_SUPPORTED_PYTHON_VERSION),
            sys_exit_code=1,
        )

    console = _set_up_logging(verbose, quiet, no_ansi)

    console.print(
        f"\n  [dark_orange]///[/] [bold][link=https://multiqc.info]MultiQC[/link][/] :mag: [dim]| v{config.version}\n"
    )
    logger.debug(f"This is MultiQC v{config.version}")

    _init(
        analysis_dir=analysis_dir,
        dirs=dirs,
        dirs_depth=dirs_depth,
        no_clean_sname=no_clean_sname,
        title=title,
        report_comment=report_comment,
        template=template,
        module=module,
        require_logs=require_logs,
        exclude=exclude,
        outdir=outdir,
        ignore=ignore,
        ignore_samples=ignore_samples,
        use_filename_as_sample_name=use_filename_as_sample_name,
        replace_names=replace_names,
        sample_names=sample_names,
        sample_filters=sample_filters,
        file_list=file_list,
        filename=filename,
        make_data_dir=make_data_dir,
        no_data_dir=no_data_dir,
        data_format=data_format,
        zip_data_dir=zip_data_dir,
        force=force,
        ignore_symlinks=ignore_symlinks,
        no_report=no_report,
        export_plots=export_plots,
        plots_flat=plots_flat,
        plots_interactive=plots_interactive,
        strict=strict,
        lint=lint,
        development=development,
        make_pdf=make_pdf,
        no_megaqc_upload=no_megaqc_upload,
        config_file=config_file,
        cl_config=cl_config,
        profile_runtime=profile_runtime,
        no_ansi=no_ansi,
        custom_css_files=custom_css_files,
        **kwargs,
    )

    try:
        run_modules, run_module_names = _file_search(
            file_list=file_list,
            ignore=ignore,
            ignore_samples=ignore_samples,
        )

        # Create the temporary working directories
        tmp_dir = tempfile.mkdtemp()
        logger.debug(f"Using temporary directory for creating report: {tmp_dir}")

        # Load the template
        template_mod = config.avail_templates[config.template].load()

        _run_modules(
            tmp_dir=tmp_dir,
            template_mod=template_mod,
            filename=filename,
            run_modules=run_modules,
            run_module_names=run_module_names,
        )

        _general_stats_table()

        _export_sources()

        _write_json_dump()

        _write_html_and_data(
            tmp_dir=tmp_dir,
            template_mod=template_mod,
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

    if report.num_mpl_plots > 0 and not config.plots_force_flat:
        if not config.plots_force_interactive:
            console.print(
                "[blue]|           multiqc[/] | "
                "Flat-image plots used. Disable with '--interactive'. "
                "See [link=https://multiqc.info/docs/#flat--interactive-plots]docs[/link]."
            )

    sys_exit_code = 0
    if config.strict and len(report.lint_errors) > 0:
        logger.error(f"Found {len(report.lint_errors)} linting errors!\n" + "\n".join(report.lint_errors))
        sys_exit_code = 1

    logger.info("MultiQC complete")

    # Move the log file into the data directory
    log.move_tmp_log(logger)

    return RunResult(sys_exit_code=sys_exit_code)


def _set_up_logging(verbose, quiet, no_ansi):
    # Set up logging level
    loglevel = log.LEVELS.get(min(verbose, 1), "INFO")
    if quiet:
        loglevel = "WARNING"
        config.quiet = True
    log.init_log(logger, loglevel=loglevel, no_ansi=no_ansi)

    console = rich.console.Console(
        stderr=True,
        highlight=False,
        force_terminal=util_functions.force_term_colors(),
        color_system=None if no_ansi else "auto",
    )

    return console


def _init(
    analysis_dir,
    dirs=False,
    dirs_depth=None,
    no_clean_sname=False,
    title=None,
    report_comment=None,
    template=None,
    module=(),
    require_logs=False,
    exclude=(),
    outdir=None,
    use_filename_as_sample_name=False,
    replace_names=None,
    sample_names=None,
    sample_filters=None,
    make_data_dir=False,
    no_data_dir=False,
    data_format=None,
    zip_data_dir=False,
    force=True,
    ignore_symlinks=False,
    no_report=False,
    export_plots=False,
    plots_flat=False,
    plots_interactive=False,
    strict=False,
    lint=False,  # Deprecated since v1.17
    development=False,
    make_pdf=False,
    no_megaqc_upload=False,
    config_file=(),
    cl_config=(),
    profile_runtime=False,
    no_ansi=False,
    custom_css_files=(),
    **kwargs,
) -> None:
    """
    Load config and set up key variables.
    """
    # Load config files
    plugin_hooks.mqc_trigger("before_config")
    config.mqc_load_userconfig(config_file)
    plugin_hooks.mqc_trigger("config_loaded")

    # Command-line config YAML
    if len(cl_config) > 0:
        config.mqc_cl_config(cl_config)

    report.init()

    # Log the command used to launch MultiQC
    report.multiqc_command = " ".join(sys.argv)
    logger.debug(f"Command used: {report.multiqc_command}")

    # Check that we're running the latest version of MultiQC
    if config.no_version_check is not True:
        try:
            # Fetch the version info from the API
            meta = {
                "version_multiqc": config.short_version,
                "version_python": platform.python_version(),
                "operating_system": platform.system(),
                "is_docker": os.path.exists("/.dockerenv"),
                "is_singularity": os.path.exists("/.singularity.d"),
                "is_conda": os.path.exists(os.path.join(sys.prefix, "conda-meta")),
                "is_ci": strtobool(os.getenv("CI", False)),
            }
            wait_seconds = 2
            try:
                r = requests.get(config.version_check_url, params=meta, timeout=wait_seconds)
            except requests.exceptions.Timeout as e:
                logger.debug(
                    f"Timed out after waiting for {wait_seconds}s for multiqc.info to check latest version: {e}"
                )
            except requests.exceptions.RequestException as e:
                logger.debug(f"Could not connect to multiqc.info for version check: {e}")
            else:
                release_info = r.json()
                # Broadcast log messages if found
                for msg in release_info.get("broadcast_messages", []):
                    if msg.get("message"):
                        level = msg.get("level")
                        if level not in ["debug", "info", "warning", "error", "critical"]:
                            level = "info"
                        getattr(logger, level)(msg["message"])
                # Available update log if newer
                remove_version = version.parse(re.sub(r"[^0-9.]", "", release_info["latest_release"]["version"]))
                this_version = version.parse(re.sub(r"[^0-9.]", "", config.short_version))
                if remove_version > this_version:
                    logger.warning(f"MultiQC Version {release_info['latest_release']['version']} now available!")
                logger.debug(
                    f"Latest MultiQC version is {release_info['latest_release']['version']}, "
                    f"released {release_info['latest_release']['release_date']}"
                )
        except Exception as e:
            logger.debug(f"Could not connect to multiqc.info for version check: {e}")

    # Set up key variables (overwrite config vars from command line)
    if template is not None:
        config.template = template
    if title is not None:
        config.title = title
    if report_comment is not None:
        config.report_comment = report_comment
    if dirs is True:
        config.prepend_dirs = dirs
    if dirs_depth is not None:
        config.prepend_dirs = True
        config.prepend_dirs_depth = dirs_depth
    # Clean up analysis_dir if a string (interactive environment only)
    if isinstance(analysis_dir, str):
        analysis_dir = [analysis_dir]
    config.analysis_dir = analysis_dir
    if outdir is not None:
        config.output_dir = os.path.realpath(outdir)
    if use_filename_as_sample_name:
        config.use_filename_as_sample_name = True
        logger.info("Using log filenames for sample names")
    if make_data_dir:
        config.make_data_dir = True
    if no_data_dir:
        config.make_data_dir = False
    if force:
        config.force = True
    if ignore_symlinks:
        config.ignore_symlinks = True
    if zip_data_dir:
        config.zip_data_dir = True
    if data_format is not None:
        config.data_format = data_format
    if export_plots:
        config.export_plots = True
    if no_report:
        config.make_report = False
    if plots_flat:
        config.plots_force_flat = True
    if plots_interactive:
        config.plots_force_interactive = True
    if lint or config.lint:  # Deprecated since v1.17
        logger.warning(
            "DEPRECIATED: The --lint option is renamed to --strict since MultiQC 1.17. "
            "The old option will be removed in future MultiQC versions, please "
            "update your command line and/or configs."
        )
        strict = True
    if strict:
        config.strict = True
        config.lint = True  # Deprecated since v1.17
        strict_helpers.run_tests()
    if development:
        config.development = True
        config.export_plots = True
        if "png" not in config.export_plot_formats:
            config.export_plot_formats.append("png")
    if make_pdf:
        config.make_pdf = True
        config.template = "simple"
    if no_megaqc_upload:
        config.megaqc_upload = False
    else:
        config.megaqc_upload = True
    if no_clean_sname:
        config.fn_clean_sample_names = False
        logger.info("Not cleaning sample names")
    if replace_names:
        config.load_replace_names(replace_names)
    if sample_names:
        config.load_sample_names(sample_names)
    config.load_show_hide(sample_filters)
    if len(module) > 0:
        config.run_modules = module
    if len(exclude) > 0:
        config.exclude_modules = exclude
    if require_logs:
        config.require_logs = True
    if profile_runtime:
        config.profile_runtime = True
    if no_ansi:
        config.no_ansi = True
    if custom_css_files:
        config.custom_css_files.extend(custom_css_files)
    config.kwargs = kwargs  # Plugin command line options

    plugin_hooks.mqc_trigger("execution_start")

    logger.debug(f"Working dir : {os.getcwd()}")
    if config.make_pdf:
        logger.info("--pdf specified. Using non-interactive HTML template.")
    logger.debug(f"Template    : {config.template}")
    if config.strict:
        logger.info(
            "Strict mode specified. Will exit early if a module or a template crashed, and will "
            "give warnings if anything is not optimally configured in a module or a template."
        )

    logger.debug("Running Python " + sys.version.replace("\n", " "))


def _file_search(
    file_list=False,
    ignore=(),
    ignore_samples=(),
) -> Tuple[List, List]:
    """
    Search log files and set up the list of modules to run.
    """

    # Add files if --file-list option is given
    if file_list:
        if len(config.analysis_dir) > 1:
            raise _RunError("If --file-list is given, analysis_dir should have only one plain text file.")
        file_list_path = config.analysis_dir[0]
        config.analysis_dir = []
        with open(file_list_path) as in_handle:
            for line in in_handle:
                if os.path.exists(line.strip()):
                    path = os.path.abspath(line.strip())
                    config.analysis_dir.append(path)
        if len(config.analysis_dir) == 0:
            raise _RunError(
                f"No files or directories were added from {file_list_path} using --file-list option."
                f"Please, check that {file_list_path} contains correct paths."
            )

    if len(ignore) > 0:
        logger.debug(f"Ignoring files, directories and paths that match: {', '.join(ignore)}")
        config.fn_ignore_files.extend(ignore)
        config.fn_ignore_dirs.extend(ignore)
        config.fn_ignore_paths.extend(ignore)
    if len(ignore_samples) > 0:
        logger.debug(f"Ignoring sample names that match: {', '.join(ignore_samples)}")
        config.sample_names_ignore.extend(ignore_samples)

    # Print some status updates
    if config.title is not None:
        logger.info(f"Report title: {config.title}")
    if config.prepend_dirs:
        logger.info("Prepending directory to sample names")

    # Prep module configs
    config.top_modules = [m if isinstance(m, dict) else {m: {}} for m in config.top_modules]
    config.module_order = [m if isinstance(m, dict) else {m: {}} for m in config.module_order]
    mod_keys = [list(m.keys())[0] for m in config.module_order]

    # Lint the module configs
    if config.strict:
        for m in config.avail_modules.keys():
            if m not in mod_keys:
                errmsg = f"LINT: Module '{m}' not found in config.module_order"
                logger.error(errmsg)
                report.lint_errors.append(errmsg)

    # Get the list of modules we want to run, in the order that we want them
    run_modules = [m for m in config.top_modules if list(m.keys())[0] in config.avail_modules.keys()]
    run_modules.extend([{m: {}} for m in config.avail_modules.keys() if m not in mod_keys and m not in run_modules])
    run_modules.extend(
        [
            m
            for m in config.module_order
            if list(m.keys())[0] in config.avail_modules.keys()
            and list(m.keys())[0] not in [list(rm.keys())[0] for rm in run_modules]
        ]
    )

    # Check the module names that were requested to run explicitly
    if len(getattr(config, "run_modules", {})) > 0:
        unknown_modules = [m for m in config.run_modules if m not in config.avail_modules.keys()]
        if unknown_modules:
            logger.error(f"Module(s) in config.run_modules are unknown: {', '.join(unknown_modules)}")
        if len(unknown_modules) == len(config.run_modules):
            raise _RunError("No available modules to run!")
        config.run_modules = [m for m in config.run_modules if m in config.avail_modules.keys()]
        run_modules = [m for m in run_modules if list(m.keys())[0] in config.run_modules]
        logger.info(f"Only using modules: {', '.join(config.run_modules)}")
    if len(getattr(config, "exclude_modules", {})) > 0:
        logger.info("Excluding modules '{}'".format("', '".join(config.exclude_modules)))
        if "general_stats" in config.exclude_modules:
            config.skip_generalstats = True
            config.exclude_modules = tuple(x for x in config.exclude_modules if x != "general_stats")
        run_modules = [m for m in run_modules if list(m.keys())[0] not in config.exclude_modules]
    if len(run_modules) == 0:
        raise _RunError("No analysis modules specified!")
    run_module_names = [list(m.keys())[0] for m in run_modules]
    logger.debug(f"Analysing modules: {', '.join(run_module_names)}")

    # Add custom content section names
    try:
        if "custom_content" in run_module_names:
            run_module_names.extend(config.custom_data.keys())
    except AttributeError:
        pass  # custom_data not in config

    # Always run software_versions module to collect version YAML files
    # Use config.skip_versions_section to exclude from report
    if "software_versions" not in run_module_names:
        run_module_names.append("software_versions")

    for d in config.analysis_dir:
        logger.info(f"Search path : {os.path.abspath(d)}")

    # FILE SEARCH. Heavy part. Go over provided paths and prepare a list of relevant log files.
    report.get_filelist(run_module_names)
    # END FILE SEARCH

    # Only run the modules for which any files were found
    non_empty_modules = {key.split("/")[0].lower() for key, files in report.files.items() if len(files) > 0}
    # Always run custom content, as it can have data purely from a MultiQC config file (no search files)
    if "custom_content" not in non_empty_modules:
        non_empty_modules.add("custom_content")
    run_modules = [m for m in run_modules if list(m.keys())[0].lower() in non_empty_modules]
    run_module_names = [list(m.keys())[0] for m in run_modules]
    if not _required_logs_found(run_module_names):
        raise _RunError()

    return run_modules, run_module_names


def _run_modules(
    tmp_dir: str,
    template_mod,
    run_modules: List[Dict[str, Optional[Dict]]],
    run_module_names: List[str],
    filename: str,
) -> None:
    config.data_tmp_dir = os.path.join(tmp_dir, "multiqc_data")
    if filename != "stdout" and config.make_data_dir is True:
        config.data_dir = config.data_tmp_dir
        os.makedirs(config.data_dir)
    else:
        config.data_dir = None
    config.plots_tmp_dir = os.path.join(tmp_dir, "multiqc_plots")
    if filename != "stdout" and config.export_plots is True:
        config.plots_dir = config.plots_tmp_dir
        os.makedirs(config.plots_dir)
    else:
        config.plots_dir = None
    if not config.make_report:
        config.output_fn = None

    # Add an output subdirectory if specified by template
    try:
        config.output_dir = os.path.join(config.output_dir, template_mod.output_subdir)
    except AttributeError:
        pass  # No subdirectory variable given

    # Run the modules!
    plugin_hooks.mqc_trigger("before_modules")
    report.modules_output = list()
    sys_exit_code = 0
    total_mods_starttime = time.time()
    for mod_idx, mod_dict in enumerate(run_modules):
        mod_starttime = time.time()
        this_module: str = list(mod_dict.keys())[0]
        mod_cust_config: Dict = list(mod_dict.values())[0] or {}
        # noinspection PyBroadException
        try:
            entry_point: EntryPoint = config.avail_modules[this_module]
            module_initializer: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
            module_initializer.mod_cust_config = mod_cust_config

            # *********************************************
            # RUN MODULE. Heavy part. Run module logic to parse logs and prepare plot data.
            modules = module_initializer()
            # END RUN MODULE
            # *********************************************

            if not isinstance(modules, list):
                modules = [modules]

            report.modules_output.extend(modules)

            if config.make_report:
                # Copy over css & js files if requested by the theme
                try:
                    for to, path in report.modules_output[-1].css.items():
                        copy_to = os.path.join(tmp_dir, to)
                        os.makedirs(os.path.dirname(copy_to))
                        shutil.copyfile(path, copy_to)
                except OSError as e:
                    if e.errno == errno.EEXIST:
                        pass
                    else:
                        raise
                except AttributeError:
                    pass
                try:
                    for to, path in report.modules_output[-1].js.items():
                        copy_to = os.path.join(tmp_dir, to)
                        os.makedirs(os.path.dirname(copy_to))
                        shutil.copyfile(path, copy_to)
                except OSError as e:
                    if e.errno == errno.EEXIST:
                        pass
                    else:
                        raise
                except AttributeError:
                    pass

        except ModuleNoSamplesFound:
            logger.debug(f"No samples found: {this_module}")
        except UserWarning:  # UserWarning deprecated from 1.16
            msg = f"DEPRECIATED: Please raise 'ModuleNoSamplesFound' instead of 'UserWarning' in module: {this_module}"
            if config.strict:
                logger.error(msg)
                report.lint_errors.append(msg)
            else:
                logger.debug(msg)
            logger.debug(f"No samples found: {this_module}")
        except KeyboardInterrupt:
            shutil.rmtree(tmp_dir)
            logger.critical(
                "User Cancelled Execution!\n{eq}\n{tb}{eq}\n".format(eq=("=" * 60), tb=traceback.format_exc())
                + "User Cancelled Execution!\nExiting MultiQC..."
            )
            sys.exit(1)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive debugging of modules.
                raise

            # Flag the error, but carry on
            class CustomTraceback:
                def __rich_console__(self, console: rich.console.Console, options: rich.console.ConsoleOptions):
                    sys_tb = sys.exc_info()
                    issue_url = "https://github.com/MultiQC/MultiQC/issues/new?template=bug_report.md&title={}%20module%20-%20{}".format(
                        this_module, sys_tb[0].__name__
                    )
                    yield (
                        "Please copy this log and report it at [bright_blue][link={}]https://github.com/MultiQC/MultiQC/issues[/link][/] \n"
                        "[bold underline]Please attach a file that triggers the error.[/] The last file found was: [green]{}[/]\n".format(
                            issue_url, report.last_found_file
                        )
                    )
                    yield Syntax(traceback.format_exc(), "python")

                def __rich_measure__(self, console: rich.console.Console, options: rich.console.ConsoleOptions):
                    tb_width = max([len(line) for line in traceback.format_exc().split("\n")])
                    try:
                        log_width = 71 + len(report.last_found_file)
                    except TypeError:
                        log_width = 71
                    panel_width = max(tb_width, log_width)
                    return rich.console.Measurement(panel_width, panel_width)

            console = rich.console.Console(
                stderr=True,
                force_terminal=util_functions.force_term_colors(),
                color_system=None if config.no_ansi else "auto",
            )
            console.print(
                rich.panel.Panel(
                    CustomTraceback(),
                    title=f"Oops! The '[underline]{this_module}[/]' MultiQC module broke...",
                    expand=False,
                    border_style="red",
                    style="on #272822",
                )
            )
            # Still log.debug this so that it ends up in the log file - above is just stderr for now
            logger.debug(
                f"Oops! The '{this_module}' MultiQC module broke...\n"
                + ("=" * 80)
                + "\n"
                + traceback.format_exc()
                + ("=" * 80)
            )
            # Exit code 1 for CI failures etc
            sys_exit_code = 1

        report.runtimes["mods"][run_module_names[mod_idx]] = time.time() - mod_starttime
    report.runtimes["total_mods"] = time.time() - total_mods_starttime

    # Again, if config.require_logs is set, check if for all explicitly requested
    # modules samples were found.
    if not _required_logs_found([m.anchor for m in report.modules_output]):
        raise _RunError()

    # Update report with software versions provided in configs
    software_versions.update_versions_from_config(config, report)

    # Add section for software versions if any are found
    if not config.skip_versions_section and report.software_versions:
        # Importing here to avoid circular imports
        from multiqc.modules.software_versions import MultiqcModule

        report.modules_output.append(MultiqcModule())

    # Special-case module if we want to profile the MultiQC running time
    if config.profile_runtime:
        from multiqc.modules.profile_runtime import MultiqcModule

        report.modules_output.append(MultiqcModule())

    # Did we find anything?
    if len(report.modules_output) == 0:
        logger.warning("No analysis results found. Cleaning up…")
        shutil.rmtree(tmp_dir)
        logger.info("MultiQC complete")
        # Exit with an error code if a module broke
        raise _RunError(sys_exit_code=sys_exit_code)

    if config.make_report:
        # Sort the report module output if we have a config
        if len(getattr(config, "report_section_order", {})) > 0:
            section_id_order = {}
            idx = 10
            for mod in reversed(report.modules_output):
                section_id_order[mod.anchor] = idx
                idx += 10
            for anchor, ss in config.report_section_order.items():
                if anchor not in section_id_order.keys():
                    logger.debug(f"Reordering sections: anchor '{anchor}' not found.")
                    continue
                if ss.get("order") is not None:
                    section_id_order[anchor] = ss["order"]
                if ss.get("after") in section_id_order.keys():
                    section_id_order[anchor] = section_id_order[ss["after"]] + 1
                if ss.get("before") in section_id_order.keys():
                    section_id_order[anchor] = section_id_order[ss["before"]] - 1
            sorted_ids = sorted(section_id_order, key=section_id_order.get)
            report.modules_output = [
                mod for i in reversed(sorted_ids) for mod in report.modules_output if mod.anchor == i
            ]

        # Sort the report sections if we have a config
        # Basically the same as above, but sections within a module
        if len(getattr(config, "report_section_order", {})) > 0:
            # Go through each module
            for midx, mod in enumerate(report.modules_output):
                section_id_order = {}
                # Get a list of the section anchors
                idx = 10
                for s in mod.sections:
                    section_id_order[s["anchor"]] = idx
                    idx += 10
                # Go through each section to be reordered
                for anchor, ss in config.report_section_order.items():
                    # Section to be moved is not in this module
                    if anchor not in section_id_order.keys():
                        logger.debug(f"Reordering sections: anchor '{anchor}' not found for module '{mod.name}'.")
                        continue
                    if ss == "remove":
                        section_id_order[anchor] = False
                        continue
                    if ss.get("order") is not None:
                        section_id_order[anchor] = ss["order"]
                    if ss.get("after") in section_id_order.keys():
                        section_id_order[anchor] = section_id_order[ss["after"]] + 1
                    if ss.get("before") in section_id_order.keys():
                        section_id_order[anchor] = section_id_order[ss["before"]] - 1
                # Remove module sections
                section_id_order = {s: o for s, o in section_id_order.items() if o is not False}
                # Sort the module sections
                sorted_ids = sorted(section_id_order, key=section_id_order.get)
                report.modules_output[midx].sections = [s for i in sorted_ids for s in mod.sections if s["anchor"] == i]

    plugin_hooks.mqc_trigger("after_modules")


def _general_stats_table() -> None:
    """
    Construct HTML for the general stats table.
    """

    # Remove empty data sections from the General Stats table
    empty_keys = [i for i, d in enumerate(report.general_stats_data[:]) if len(d) == 0]
    empty_keys.sort(reverse=True)
    for i in empty_keys:
        del report.general_stats_data[i]
        del report.general_stats_headers[i]

    # Add general-stats IDs to table row headers
    for idx, h in enumerate(report.general_stats_headers):
        for k in h.keys():
            if "rid" not in h[k]:
                h[k]["rid"] = re.sub(r"\W+", "_", k).strip().strip("_")
            ns_html = re.sub(r"\W+", "_", h[k]["namespace"]).strip().strip("_").lower()
            report.general_stats_headers[idx][k]["rid"] = report.save_htmlid(
                f"mqc-generalstats-{ns_html}-{h[k]['rid']}"
            )

    all_hidden = True
    for headers in report.general_stats_headers:
        for h in headers.values():
            if not h.get("hidden", False):
                all_hidden = False
                break

    # Generate the General Statistics HTML & write to file
    if len(report.general_stats_data) > 0 and not config.skip_generalstats and not all_hidden:
        pconfig = {
            "id": "general_stats_table",
            "table_title": "General Statistics",
            "save_file": True,
            "raw_data_fn": "multiqc_general_stats",
        }
        report.general_stats_html = table.plot(report.general_stats_data, report.general_stats_headers, pconfig)
    else:
        config.skip_generalstats = True


def _export_sources() -> None:
    """
    Dump data sources.
    """

    if config.data_dir is not None:
        # Write the report sources to disk
        report.data_sources_tofile()

        # Create a file with the module DOIs
        report.dois_tofile(report.modules_output)


def _write_json_dump() -> None:
    """
    Write JSON with plot data. Useful for MegaQC and for loading in interactive environments.
    """

    plugin_hooks.mqc_trigger("before_report_generation")

    # Data Export / MegaQC integration - save report data to file or send report data to an API endpoint
    if config.data_dump_file or (config.megaqc_url and config.megaqc_upload):
        multiqc_json_dump = util_functions.multiqc_dump_json(report)
        if config.data_dump_file:
            util_functions.write_data_file(multiqc_json_dump, "multiqc_data", False, "json")
        if config.megaqc_url:
            megaqc.multiqc_api_post(multiqc_json_dump)

    if config.development:
        with open(os.path.join(config.data_dir, "multiqc_plots.js"), "w") as f:
            f.write(json.dumps(report.plot_data))


def _write_html_and_data(
    tmp_dir: str,
    template_mod,
    filename: str,
) -> None:
    """
    Make the final report path & data directories
    """

    if filename == "stdout":
        config.output_fn = sys.stdout
        logger.info("Printing report to stdout")
    else:
        if filename is not None and filename.endswith(".html"):
            filename = filename[:-5]
        if filename is None and config.title is not None:
            filename = re.sub(r"[^\w.-]", "", re.sub(r"[-\s]+", "-", config.title)).strip()
            filename += "_multiqc_report"
        if filename is not None:
            if "output_fn_name" not in config.nondefault_config:
                config.output_fn_name = f"{filename}.html"
            if "data_dir_name" not in config.nondefault_config:
                config.data_dir_name = f"{filename}_data"
            if "plots_dir_name" not in config.nondefault_config:
                config.plots_dir_name = f"{filename}_plots"
        if not config.output_fn_name.endswith(".html"):
            config.output_fn_name = f"{config.output_fn_name}.html"

        if config.make_report:
            config.output_fn = os.path.join(config.output_dir, config.output_fn_name)
        config.data_dir = os.path.join(config.output_dir, config.data_dir_name)
        config.plots_dir = os.path.join(config.output_dir, config.plots_dir_name)
        # del config.data_dir_name
        # del config.plots_dir_name
        deleted_report = False
        deleted_data_dir = False
        deleted_export_plots = False
        # Check for existing reports and remove if -f was specified
        if (
            (config.make_report and os.path.exists(config.output_fn))
            or (config.make_data_dir and os.path.exists(config.data_dir))
            or (config.export_plots and os.path.exists(config.plots_dir))
        ):
            if config.force:
                if config.make_report and os.path.exists(config.output_fn):
                    deleted_report = True
                    os.remove(config.output_fn)
                if config.make_data_dir and os.path.exists(config.data_dir):
                    deleted_data_dir = True
                    shutil.rmtree(config.data_dir)
                if config.export_plots and os.path.exists(config.plots_dir):
                    deleted_export_plots = True
                    shutil.rmtree(config.plots_dir)
            else:
                # Set up the base names of the report and the data dir
                report_num = 1
                dir_base = os.path.basename(config.data_dir)
                plots_base = os.path.basename(config.plots_dir)

                # Iterate through appended numbers until we find one that's free
                while (
                    (config.make_report and os.path.exists(config.output_fn))
                    or (config.make_data_dir and os.path.exists(config.data_dir))
                    or (config.export_plots and os.path.exists(config.plots_dir))
                ):
                    if config.make_report:
                        report_base, report_ext = os.path.splitext(config.output_fn_name)
                        config.output_fn = os.path.join(config.output_dir, f"{report_base}_{report_num}{report_ext}")
                    config.data_dir = os.path.join(config.output_dir, f"{dir_base}_{report_num}")
                    config.plots_dir = os.path.join(config.output_dir, f"{plots_base}_{report_num}")
                    report_num += 1
                if config.make_report:
                    config.output_fn_name = os.path.basename(config.output_fn)
                config.data_dir_name = os.path.basename(config.data_dir)
                config.plots_dir_name = os.path.basename(config.plots_dir)
                logger.info("Existing reports found, adding suffix to filenames. Use '--force' to overwrite.")

        # Make directories for report if needed
        if config.make_report:
            if not os.path.exists(os.path.dirname(config.output_fn)):
                os.makedirs(os.path.dirname(config.output_fn))
            logger.info(
                "Report      : {}{}".format(
                    os.path.relpath(config.output_fn),
                    "   (overwritten)" if deleted_report else "",
                )
            )
        else:
            logger.info("Report      : None")

        if config.make_data_dir is False:
            logger.info("Data        : None")
        else:
            # Make directories for data_dir
            logger.info(
                "Data        : {}{}".format(
                    os.path.relpath(config.data_dir),
                    "   (overwritten)" if deleted_data_dir else "",
                )
            )
            # Modules have run, so data directory should be complete by now. Move its contents.
            logger.debug(f"Moving data file from '{config.data_tmp_dir}' to '{config.data_dir}'")
            shutil.copytree(
                config.data_tmp_dir,
                config.data_dir,
                # Override default shutil.copy2 function to copy files. The default
                # function copies times and mode, which we want to avoid on purpose
                # to get around the problem with mounted CIFS shares (see #625).
                # shutil.copyfile only copies the file without any metadata.
                copy_function=shutil.copyfile,
            )
            shutil.rmtree(config.data_tmp_dir)

        if config.output_fn is not None:
            logger.debug(f"Full report path: {os.path.realpath(config.output_fn)}")

        # Copy across the static plot images if requested
        if config.export_plots:
            config.plots_dir = os.path.join(config.output_dir, config.plots_dir_name)
            if os.path.exists(config.plots_dir):
                if config.force:
                    shutil.rmtree(config.plots_dir)
                else:
                    logger.error(f"Output directory {config.plots_dir} already exists.")
                    logger.info("Use -f or --force to overwrite existing reports")
                    shutil.rmtree(tmp_dir)
                    raise _RunError()
            logger.info(
                "Plots       : {}{}".format(
                    os.path.relpath(config.plots_dir),
                    "   (overwritten)" if deleted_export_plots else "",
                )
            )

            # Modules have run, so plots directory should be complete by now. Move its contents.
            logger.debug(f"Moving plots directory from '{config.plots_tmp_dir}' to '{config.plots_dir}'")
            shutil.copytree(
                config.plots_tmp_dir,
                config.plots_dir,
                # Override default shutil.copy2 function to copy files. The default
                # function copies times and mode, which we want to avoid on purpose
                # to get around the problem with mounted CIFS shares (see #625).
                # shutil.copyfile only copies the file without any metadata.
                copy_function=shutil.copyfile,
            )
            shutil.rmtree(config.plots_tmp_dir)

    plugin_hooks.mqc_trigger("before_template")

    # Generate report if required
    if config.make_report:
        # Load in parent template files first if a child theme
        parent_template = None
        try:
            parent_template = config.avail_templates[template_mod.template_parent].load()
        except AttributeError:
            pass  # Not a child theme
        else:
            shutil.copytree(parent_template.template_dir, tmp_dir, dirs_exist_ok=True)

        # Copy the template files to the tmp directory (`dirs_exist_ok` makes sure
        # parent template files are overwritten)
        shutil.copytree(template_mod.template_dir, tmp_dir, dirs_exist_ok=True)

        # Function to include file contents in Jinja template
        def include_file(name, fdir=tmp_dir, b64=False):
            try:
                if fdir is None:
                    fdir = ""
                path = os.path.join(fdir, name)

                if config.development:
                    if os.path.exists(dev_path := os.path.join(template_mod.template_dir, name)):
                        fdir = template_mod.template_dir
                        name = dev_path
                        path = dev_path
                    elif parent_template and os.path.exists(
                        dev_path := os.path.join(parent_template.template_dir, name)
                    ):
                        fdir = template_mod.template_dir
                        name = dev_path
                        path = dev_path

                    if re.match(r".*\.min\.(js|css)$", name):
                        unminimized_name = re.sub(r"\.min\.", ".", name)
                        if os.path.exists(os.path.join(fdir, unminimized_name)):
                            name = unminimized_name

                    if name.endswith(".js"):
                        return f'</script><script type="text/javascript" src="{name}">'
                    if name.endswith(".css"):
                        return f'</style><link rel="stylesheet" href="{name}">'

                if b64:
                    with io.open(path, "rb") as f:
                        return base64.b64encode(f.read()).decode("utf-8")
                else:
                    with io.open(path, "r", encoding="utf-8") as f:
                        return f.read()
            except (OSError, IOError) as e:
                logger.error(f"Could not include file '{name}': {e}")

        # Load the report template
        try:
            env = jinja2.Environment(loader=jinja2.FileSystemLoader(tmp_dir))
            env.globals["include_file"] = include_file
            j_template = env.get_template(template_mod.base_fn, globals={"development": config.development})
        except:  # noqa: E722
            raise IOError(f"Could not load {config.template} template file '{template_mod.base_fn}'")

        # Compress the report plot JSON data
        runtime_compression_start = time.time()
        logger.debug("Compressing plot data")
        report.plot_compressed_json = report.compress_json(report.plot_data)
        report.runtimes["total_compression"] = time.time() - runtime_compression_start

        # Use jinja2 to render the template and overwrite
        config.analysis_dir = [os.path.realpath(d) for d in config.analysis_dir]
        report_output = j_template.render(report=report, config=config)
        if filename == "stdout":
            print(report_output.encode("utf-8"), file=sys.stdout)
        else:
            try:
                with io.open(config.output_fn, "w", encoding="utf-8") as f:
                    print(report_output, file=f)
            except IOError as e:
                raise IOError(f"Could not print report to '{config.output_fn}' - {IOError(e)}")

            # Copy over files if requested by the theme
            try:
                for f in template_mod.copy_files:
                    fn = os.path.join(tmp_dir, f)
                    dest_dir = os.path.join(os.path.dirname(config.output_fn), f)
                    shutil.copytree(fn, dest_dir, dirs_exist_ok=True)
            except AttributeError:
                pass  # No files to copy

    # Clean up temporary directory
    shutil.rmtree(tmp_dir)

    # Zip the data directory if requested
    if config.zip_data_dir and config.data_dir is not None:
        shutil.make_archive(config.data_dir, "zip", config.data_dir)
        shutil.rmtree(config.data_dir)

    # Try to create a PDF if requested
    if config.make_pdf:
        try:
            pdf_fn_name = config.output_fn.replace(".html", ".pdf")
            pandoc_call = [
                "pandoc",
                "--standalone",
                config.output_fn,
                "--output",
                pdf_fn_name,
                "--pdf-engine=xelatex",
                "-V",
                "documentclass=article",
                "-V",
                "geometry=margin=1in",
                "-V",
                "title=",
            ]
            if config.pandoc_template is not None:
                pandoc_call.append(f"--template={config.pandoc_template}")
            logger.debug(f"Attempting Pandoc conversion to PDF with following command:\n{' '.join(pandoc_call)}")
            pdf_exit_code = subprocess.call(pandoc_call)
            if pdf_exit_code != 0:
                logger.error("Error creating PDF! Pandoc returned a non-zero exit code.")
            else:
                logger.info(f"PDF Report  : {pdf_fn_name}")
        except OSError as e:
            if e.errno == errno.ENOENT:
                logger.error("Error creating PDF - pandoc not found. Is it installed? http://pandoc.org/")
            else:
                logger.error(
                    "Error creating PDF! Something went wrong when creating the PDF\n"
                    + ("=" * 60)
                    + f"\n{traceback.format_exc()}\n"
                    + ("=" * 60)
                )


def _required_logs_found(modules_with_logs):
    if config.require_logs:
        required_modules_with_no_logs = [
            m
            for m in getattr(config, "run_modules", [])
            if m.lower() not in [m.lower() for m in modules_with_logs]
            and m.lower() not in getattr(config, "exclude_modules", [])
        ]
        if required_modules_with_no_logs:
            logger.critical(
                "The following modules were explicitly requested but no log files were found: {}".format(
                    ", ".join(required_modules_with_no_logs)
                )
            )
            return False
    return True
