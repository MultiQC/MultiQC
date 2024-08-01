"""
The main function to run MultiQC.
Primarily called by multiqc.__main__.py
Imported by __init__.py so available as multiqc.run()
"""

import logging
import os
import subprocess
import sys
import time
import traceback
from typing import Tuple, Optional

import rich_click as click

from multiqc import config, report
from multiqc.core import plugin_hooks, log_and_rich
from multiqc.core.exceptions import RunError, NoAnalysisFound
from multiqc.core.exec_modules import exec_modules
from multiqc.core.file_search import file_search
from multiqc.core.order_modules_and_sections import order_modules_and_sections
from multiqc.core.update_config import update_config, ClConfig
from multiqc.core.version_check import check_version
from multiqc.core.write_results import write_results
from multiqc.validation import ConfigValidationError

logger = logging.getLogger(__name__)

start_execution_time = time.time()

NO_ANSI_FLAG = "--no-ansi"

# Configuration for rich-click CLI help
click.rich_click.COLOR_SYSTEM = "auto" if NO_ANSI_FLAG not in sys.argv else None
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
emoji = log_and_rich.choose_emoji(use_rich=True)
emoji = f" {emoji}" if emoji else ""
click.rich_click.HEADER_TEXT = (
    f"[dark_orange]///[/] [bold][link=https://multiqc.info]MultiQC[/link][/]{emoji} [dim]v{config.version}[/]"
)
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
                "--profile-memory",
                "--no-megaqc-upload",
                "--no-ansi",
                "--version",
                "--help",
            ],
        },
    ],
}


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument(
    "analysis_dir",
    type=click.Path(exists=True),
    nargs=-1,
    required=True,
    metavar="[ANALYSIS DIRECTORY]",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=None,
    help="Overwrite any existing reports",
)
@click.option(
    "-d",
    "--dirs",
    "prepend_dirs",
    is_flag=True,
    default=None,
    help="Prepend directory to sample names",
)
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
    "fn_clean_sample_names",
    flag_value=False,
    default=None,
    help="Do not clean the sample names [i](leave as full file name)[/]",
)
@click.option(
    "-i",
    "--title",
    type=str,
    help="Report title. Printed as page header, used for filename if not otherwise specified.",
)
@click.option(
    "-b",
    "--comment",
    "report_comment",
    type=str,
    help="Custom comment, will be printed at the top of the report.",
)
@click.option(
    "-n",
    "--filename",
    type=str,
    help="Report filename. Use '[yellow i]stdout[/]' to print to standard out.",
)
@click.option(
    "-o",
    "--outdir",
    "output_dir",
    type=str,
    help="Create report in the specified output directory.",
)
@click.option(
    "-t",
    "--template",
    type=click.Choice(list(config.avail_templates.keys())),
    metavar=None,
    help="Report template to use.",
)
@click.option(
    "-x",
    "--ignore",
    type=str,
    multiple=True,
    metavar="GLOB EXPRESSION",
    help="Ignore analysis files",
)
@click.option(
    "--ignore-samples",
    "ignore_samples",
    type=str,
    multiple=True,
    metavar="GLOB EXPRESSION",
    help="Ignore sample names",
)
@click.option(
    "--ignore-symlinks",
    "ignore_symlinks",
    is_flag=True,
    default=None,
    help="Ignore symlinked directories and files",
)
@click.option(
    "--fn_as_s_name",
    "use_filename_as_sample_name",
    is_flag=True,
    default=None,
    help="Use the log filename as the sample name",
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
    "-l",
    "--file-list",
    "file_list",
    is_flag=True,
    default=None,
    help="Supply a file containing a list of file paths to be searched, one per row",
)
@click.option(
    "-e",
    "--exclude",
    "exclude_modules",
    metavar="[MODULE NAME]",
    type=click.Choice(sorted(["general_stats"] + list(config.avail_modules.keys()))),
    multiple=True,
    help="Do not use this module. Can specify multiple times.",
)
@click.option(
    "-m",
    "--module",
    "run_modules",
    metavar="[MODULE NAME]",
    type=click.Choice(sorted(config.avail_modules.keys())),
    multiple=True,
    help="Use only this module. Can specify multiple times.",
)
@click.option(
    "--require-logs",
    "require_logs",
    is_flag=True,
    default=None,
    help="Require all explicitly requested modules to have log files. If not, MultiQC will exit with an error.",
)
@click.option(
    "--data-dir/--no-data-dir",
    "make_data_dir",
    is_flag=True,
    default=None,
    help="Force the parsed data directory to be created.",
)
@click.option(
    "-k",
    "--data-format",
    "data_format",
    type=click.Choice(list(config.data_format_extensions.keys())),
    help="Output parsed data in a different format.",
)
@click.option(
    "-z",
    "--zip-data-dir",
    "zip_data_dir",
    is_flag=True,
    default=None,
    help="Compress the data directory.",
)
@click.option(
    "--no-report",
    "make_report",
    flag_value=False,
    default=None,
    help="Do not generate a report, only export data and plots",
)
@click.option(
    "-p",
    "--export",
    "export_plots",
    is_flag=True,
    default=None,
    help="Export plots as static images in addition to the report",
)
@click.option(
    "-fp",
    "--flat",
    "plots_force_flat",
    is_flag=True,
    default=None,
    help="Use only flat plots [i](static images)[/]",
)
@click.option(
    "-ip",
    "--interactive",
    "plots_force_interactive",
    is_flag=True,
    default=None,
    help="Use only interactive plots [i](in-browser Javascript)[/]",
)
@click.option(
    "--strict",
    "strict",
    is_flag=True,
    default=None,
    help="Don't catch exceptions, run additional code checks to help development.",
)
@click.option(
    "--development",
    "--dev",
    "development",
    is_flag=True,
    default=None,
    help="Development mode. Do not compress and minimise JS, export uncompressed plot data",
)
@click.option(
    "--pdf",
    "make_pdf",
    is_flag=True,
    default=None,
    help="Creates PDF report with the [i]'simple'[/] template. Requires [link=https://pandoc.org/]Pandoc[/] to be installed.",
)
@click.option(
    "--no-megaqc-upload",
    "no_megaqc_upload",
    is_flag=True,
    default=None,
    help="Don't upload generated report to MegaQC, even if MegaQC options are found",
)
@click.option(
    "-c",
    "--config",
    "config_files",
    type=click.Path(exists=True, readable=True),
    multiple=True,
    help="Specific config file to load, after those in MultiQC dir / home dir / working dir.",
)
@click.option(
    "--cl-config",
    type=str,
    multiple=True,
    help="Specify MultiQC config YAML on the command line",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=0,
    help="Increase output verbosity.",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    default=None,
    help="Only show log warnings",
)
@click.option(
    "--profile-runtime",
    "profile_runtime",
    is_flag=True,
    default=None,
    help="Add analysis of how long MultiQC takes to run to the report",
)
@click.option(
    "--profile-memory",
    "profile_memory",
    is_flag=True,
    default=None,
    help="Add analysis of how much memory each module uses. Note that tracking memory will increase the runtime, so the runtime metrics could scale up a few times",
)
@click.option(
    NO_ANSI_FLAG,
    is_flag=True,
    default=None,
    help="Disable coloured log output",
)
@click.option(
    "--custom-css-file",
    "custom_css_files",
    type=click.Path(exists=True, readable=True),
    multiple=True,
    help="Custom CSS file to add to the final report",
)
@click.option(
    "--clean-up/--no-clean-up",
    "clean_up",
    is_flag=True,
    default=True,
    help="Remove the temporary directory and log file after finishing",
)
@click.option(
    "--no-version-check",
    "no_version_check",
    is_flag=True,
    default=None,
    help="Disable checking the latest MultiQC version on the server",
)
@click.version_option(config.version, prog_name="multiqc")
def run_cli(analysis_dir: Tuple[str], clean_up: bool, **kwargs):
    # Main MultiQC run command for use with the click command line, complete with all click function decorators.
    # To make it easy to use MultiQC within notebooks and other locations that don't need click, we simply pass the
    # parsed variables on to a vanilla python function.

    """MultiQC aggregates results from bioinformatics analyses across many samples into a single report.

    It searches a given directory for analysis logs and compiles an HTML report.
    It's a general use tool, perfect for summarising the output from numerous
    bioinformatics tools.

    To run, supply with one or more directory to scan for analysis results.
    For example, to run in the current working directory, use '[blue bold]multiqc .[/]'
    """

    cl_config_kwargs = {k: v for k, v in kwargs.items() if k in ClConfig.model_fields}
    other_fields = {k: v for k, v in kwargs.items() if k not in ClConfig.model_fields}
    cfg = ClConfig(**cl_config_kwargs, unknown_options=other_fields)

    # Pass on to a regular function that can be used easily without click
    result = run(*analysis_dir, clean_up=clean_up, cfg=cfg, interactive=False)

    # End execution using the exit code returned from MultiQC
    sys.exit(result.sys_exit_code)


class RunResult:
    """
    Returned by a MultiQC run for interactive use. Contains the following information:

    * appropriate error code (e.g. 1 if a module broke, 0 on success)
    * error message if a module broke
    """

    def __init__(self, sys_exit_code: int = 0, message: str = ""):
        self.sys_exit_code = sys_exit_code
        self.message = message


def run(*analysis_dir, clean_up: bool = True, cfg: Optional[ClConfig] = None, interactive: bool = True) -> RunResult:
    """
    MultiQC aggregates results from bioinformatics analyses across many samples into a single report.

    It searches a given directory for analysis logs and compiles an HTML report.
    It's a general use tool, perfect for summarising the output from numerous
    bioinformatics tools.

    To run, supply with one or more directory to scan for analysis results.
    To run here, use 'multiqc .'

    See http://multiqc.info for more details.
    """

    # In case if run() is called multiple times in the same session:
    report.reset()
    config.reset()
    update_config(*analysis_dir, cfg=cfg, log_to_file=not interactive)

    check_version()

    logger.debug(f"Working dir : {os.getcwd()}")
    if config.make_pdf:
        _check_pdf_export_possible()

    logger.debug(f"Template    : {config.template}")
    if config.strict:
        logger.info(
            "Strict mode specified. Will exit early if a module or a template crashed, and will "
            "give warnings if anything is not optimally configured in a module or a template."
        )

    report.multiqc_command = " ".join(sys.argv)
    logger.debug(f"Command used: {report.multiqc_command}")

    try:
        mod_dicts_in_order = file_search()

        exec_modules(mod_dicts_in_order)

        order_modules_and_sections()

        write_results()

    except NoAnalysisFound as e:
        logger.warning(f"{e.message}. Cleaning up…")
        return RunResult(message="No analysis results found", sys_exit_code=e.sys_exit_code)

    except ConfigValidationError as e:
        logger.warning("Config validation error. Exiting because strict mode is enabled. Cleaning up…")
        return RunResult(message=e.message, sys_exit_code=1)

    except RunError as e:
        if e.message:
            logger.critical(e.message)
        return RunResult(message=e.message, sys_exit_code=e.sys_exit_code)

    except KeyboardInterrupt:
        logger.critical(
            "User Cancelled Execution!\n{eq}\n{tb}{eq}\n".format(eq=("=" * 60), tb=traceback.format_exc())
            + "User Cancelled Execution!\nExiting MultiQC..."
        )
        return RunResult(sys_exit_code=1)

    else:
        plugin_hooks.mqc_trigger("execution_finish")

        report.runtimes.total = time.time() - start_execution_time
        if config.profile_runtime:
            logger.warning(f"Run took {report.runtimes.total:.2f} seconds")
            logger.warning(f" - {report.runtimes.total_sp:.2f}s: Searching files")
            logger.warning(f" - {report.runtimes.total_mods:.2f}s: Running modules")
            if config.make_report:
                logger.warning(f" - {report.runtimes.total_compression:.2f}s: Compressing report data")
                logger.info("For more information, see the 'Run Time' section in the report")

        if report.num_flat_plots > 0 and not config.plots_force_flat:
            if not config.plots_force_interactive:
                log_and_rich.rich_console_print(
                    "[blue]|           multiqc[/] | "
                    "Flat-image plots used. Disable with '--interactive'. "
                    "See [link=https://multiqc.info/docs/#flat--interactive-plots]docs[/link]."
                )

        sys_exit_code = 0
        if config.strict and len(report.lint_errors) > 0:
            logger.error(f"Found {len(report.lint_errors)} linting errors!\n" + "\n".join(report.lint_errors))
            sys_exit_code = 1

        logger.info("MultiQC complete")
        return RunResult(sys_exit_code=sys_exit_code)

    finally:
        if clean_up:
            report.remove_tmp_dir()


def _check_pdf_export_possible():
    if subprocess.call(["which", "pandoc"]) != 0:
        logger.error(
            "`pandoc` and `pdflatex` tools are required to create a PDF report. Please install those and try "
            "again. See http://pandoc.org/installing.html for the `pandoc` installation instructions "
            "(e.g. `brew install pandoc` on macOS), and install LaTeX for `pdflatex` (e.g. `brew install basictex`"
            "on macOS). Alternatively, omit the `--pdf` option or unset `make_pdf: true` in the MultiQC config."
        )
        return RunResult(message="Pandoc is required to create PDF reports", sys_exit_code=1)

    if subprocess.call(["which", "pdflatex"]) != 0:
        logger.error(
            "The `pdflatex` tool is required to create a PDF report. Please install LaTeX and try again, "
            "e.g. `brew install basictex` on macOS. Alternatively, omit the `--pdf` option"
            "or unset `make_pdf: true` in the MultiQC config."
        )
        return RunResult(message="LaTeX is required to create PDF reports", sys_exit_code=1)

    logger.info("--pdf specified. Using non-interactive HTML template.")
