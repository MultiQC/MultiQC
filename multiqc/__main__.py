#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
multiqc.__main__
~~~~~~~~~~~~~~~~~~~~~
The main function to run MultiQC.
python multiqc.__main__.py

Also available as:
python -m multiqc .
"""
import sys

import pkg_resources
import rich_click as click

from multiqc import run

from .utils import config, util_functions

# Set up logging
logger = config.logger

# Configuration for rich-click CLI help
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.HEADER_TEXT = (
    f"[blue]/[/][green]/[/][red]/[/] [bold][link=https://multiqc.info]MultiQC[/link][/] :mag: [dim]| v{config.version}"
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
                "--tag",
                "--view-tags",
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
                "--lint",
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
@click.option("--tag", "module_tag", type=str, multiple=True, help="Use only modules which tagged with this keyword")
@click.option(
    "--view-tags",
    is_flag=True,
    callback=util_functions.view_all_tags,
    expose_value=False,
    is_eager=True,
    help="View the available tags and which modules they load",
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
@click.option("--data-dir", "make_data_dir", is_flag=True, help="Force the parsed data directory to be created.")
@click.option(
    "--no-data-dir", "no_data_dir", is_flag=True, help="Prevent the parsed data directory from being created."
)
@click.option(
    "-k",
    "--data-format",
    "data_format",
    type=click.Choice(config.data_format_extensions.keys()),
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
@click.option("--lint", "lint", is_flag=True, help="Use strict linting (validation) to help code development")
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

    # Add any extra plugin command line options
    for entry_point in pkg_resources.iter_entry_points("multiqc.cli_options.v1"):
        opt_func = entry_point.load()
        multiqc_run = opt_func(multiqc_run)

    # End execution using the exit code returned from MultiQC
    sys.exit(multiqc_run["sys_exit_code"])


if __name__ == "__main__":
    run_cli(prog_name="multiqc")
