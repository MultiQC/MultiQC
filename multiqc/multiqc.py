#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
multiqc.multiqc
~~~~~~~~~~~~~~~~~~~~~
The main function to run MultiQC. Sorry about the messy namespace.
Primarily called by multiqc.__main__.py
Imported by __init__.py so available as multiqc.run()
"""
from __future__ import print_function

import base64
import click
from distutils import version
from distutils.dir_util import copy_tree
import errno
import io
import jinja2
import os
import re
import rich
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

from rich.syntax import Syntax

try:
    # Python 3 imports
    from urllib.request import urlopen
except ImportError:
    # Python 2 imports
    from urllib2 import urlopen

    # Use UTF-8 encoding by default
    reload(sys)
    sys.setdefaultencoding("utf8")

from .plots import table
from .utils import report, plugin_hooks, megaqc, util_functions, lint_helpers, config, log, mqc_colour

start_execution_time = time.time()
logger = config.logger


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument("analysis_dir", type=click.Path(exists=True), nargs=-1, required=True, metavar="<analysis directory>")
@click.option("-f", "--force", is_flag=True, help="Overwrite any existing reports")
@click.option("-d", "--dirs", is_flag=True, help="Prepend directory to sample names")
@click.option(
    "-dd",
    "--dirs-depth",
    "dirs_depth",
    type=int,
    help="Prepend [INT] directories to sample names. Negative number to take from start of path.",
)
@click.option(
    "-s", "--fullnames", "no_clean_sname", is_flag=True, help="Do not clean the sample names (leave as full file name)"
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
@click.option("-n", "--filename", type=str, help="Report filename. Use 'stdout' to print to standard out.")
@click.option("-o", "--outdir", type=str, help="Create report in the specified output directory.")
@click.option("-t", "--template", type=click.Choice(config.avail_templates), help="Report template to use.")
@click.option(
    "--tag", "module_tag", type=str, multiple=True, help="Use only modules which tagged with this keyword, eg. RNA"
)
@click.option(
    "--view-tags",
    "--view_tags",
    is_flag=True,
    callback=util_functions.view_all_tags,
    expose_value=False,
    is_eager=True,
    help="View the available tags and which modules they load",
)
@click.option("-x", "--ignore", type=str, multiple=True, help="Ignore analysis files (glob expression)")
@click.option(
    "--ignore-samples", "ignore_samples", type=str, multiple=True, help="Ignore sample names (glob expression)"
)
@click.option("--ignore-symlinks", "ignore_symlinks", is_flag=True, help="Ignore symlinked directories and files")
@click.option(
    "--fn_as_s_name", "use_filename_as_sample_name", is_flag=True, help="Use the log filename as the sample name"
)
@click.option(
    "--replace-names",
    "replace_names",
    type=click.Path(exists=True, readable=True),
    help="Path to TSV file to rename sample names during report generation",
)
@click.option(
    "--sample-names",
    "sample_names",
    type=click.Path(exists=True, readable=True),
    help="Path to TSV file containing alternative sample names for renaming buttons in the report",
)
@click.option(
    "--sample-filters",
    "sample_filters",
    type=click.Path(exists=True, readable=True),
    help="Path to TSV file containing show/hide patterns for the report",
)
@click.option(
    "-l", "--file-list", is_flag=True, help="Supply a file containing a list of file paths to be searched, one per row"
)
@click.option(
    "-e",
    "--exclude",
    metavar="[module name]",
    type=click.Choice(sorted(["general_stats"] + list(config.avail_modules.keys()))),
    multiple=True,
    help="Do not use this module. Can specify multiple times.",
)
@click.option(
    "-m",
    "--module",
    metavar="[module name]",
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
    help="Output parsed data in a different format. Default: {}".format(config.data_format),
)
@click.option("-z", "--zip-data-dir", "zip_data_dir", is_flag=True, help="Compress the data directory.")
@click.option(
    "-p", "--export", "export_plots", is_flag=True, help="Export plots as static images in addition to the report"
)
@click.option("-fp", "--flat", "plots_flat", is_flag=True, help="Use only flat plots (static images)")
@click.option(
    "-ip", "--interactive", "plots_interactive", is_flag=True, help="Use only interactive plots (HighCharts Javascript)"
)
@click.option("--lint", "lint", is_flag=True, help="Use strict linting (validation) to help code development")
@click.option(
    "--pdf",
    "make_pdf",
    is_flag=True,
    help="Creates PDF report with 'simple' template. Requires Pandoc to be installed.",
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
@click.option(
    "--cl-config", "--cl_config", type=str, multiple=True, help="Specify MultiQC config YAML on the command line"
)
@click.option("-v", "--verbose", count=True, default=0, help="Increase output verbosity.")
@click.option("-q", "--quiet", is_flag=True, help="Only show log warnings")
@click.option("--profile-runtime", is_flag=True, help="Add analysis of how long MultiQC takes to run to the report")
@click.option("--no-ansi", is_flag=True, help="Disable coloured log output")
@click.version_option(config.version, prog_name="multiqc")
def run_cli(
    analysis_dir,
    dirs,
    dirs_depth,
    no_clean_sname,
    title,
    report_comment,
    template,
    module_tag,
    module,
    exclude,
    outdir,
    ignore,
    ignore_samples,
    use_filename_as_sample_name,
    replace_names,
    sample_names,
    sample_filters,
    file_list,
    filename,
    make_data_dir,
    no_data_dir,
    data_format,
    zip_data_dir,
    force,
    ignore_symlinks,
    export_plots,
    plots_flat,
    plots_interactive,
    lint,
    make_pdf,
    no_megaqc_upload,
    config_file,
    cl_config,
    verbose,
    quiet,
    profile_runtime,
    no_ansi,
    **kwargs,
):
    # Main MultiQC run command for use with the click command line, complete with all click function decorators.
    # To make it easy to use MultiQC within notebooks and other locations that don't need click, we simply pass the
    # parsed variables on to a vanilla python function.

    """MultiQC aggregates results from bioinformatics analyses across many samples into a single report.

    It searches a given directory for analysis logs and compiles a HTML report.
    It's a general use tool, perfect for summarising the output from numerous
    bioinformatics tools.

    To run, supply with one or more directory to scan for analysis results.
    To run here, use 'multiqc .'

    See http://multiqc.info for more details.

    Author: Phil Ewels (http://phil.ewels.co.uk)
    """

    # Use keyword arguments in case things get rearranged in the future
    multiqc_run = run(
        analysis_dir=analysis_dir,
        dirs=dirs,
        dirs_depth=dirs_depth,
        no_clean_sname=no_clean_sname,
        title=title,
        report_comment=report_comment,
        template=template,
        module_tag=module_tag,
        module=module,
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
        export_plots=export_plots,
        plots_flat=plots_flat,
        plots_interactive=plots_interactive,
        lint=lint,
        make_pdf=make_pdf,
        no_megaqc_upload=no_megaqc_upload,
        config_file=config_file,
        cl_config=cl_config,
        verbose=verbose,
        quiet=quiet,
        profile_runtime=profile_runtime,
        no_ansi=no_ansi,
        kwargs=kwargs,
    )

    # End execution using the exit code returned from MultiQC
    sys.exit(multiqc_run["sys_exit_code"])


# Main function that runs MultQC. Available to use within an interactive Python environment
def run(
    analysis_dir,
    dirs=False,
    dirs_depth=None,
    no_clean_sname=False,
    title=None,
    report_comment=None,
    template=None,
    module_tag=(),
    module=(),
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
    export_plots=False,
    plots_flat=False,
    plots_interactive=False,
    lint=False,
    make_pdf=False,
    no_megaqc_upload=False,
    config_file=(),
    cl_config=(),
    verbose=0,
    quiet=False,
    profile_runtime=False,
    no_ansi=False,
    kwargs={},
):
    """MultiQC aggregates results from bioinformatics analyses across many samples into a single report.

    It searches a given directory for analysis logs and compiles a HTML report.
    It's a general use tool, perfect for summarising the output from numerous
    bioinformatics tools.

    To run, supply with one or more directory to scan for analysis results.
    To run here, use 'multiqc .'

    See http://multiqc.info for more details.

    Author: Phil Ewels (http://phil.ewels.co.uk)
    """

    # Set up logging level
    loglevel = log.LEVELS.get(min(verbose, 1), "INFO")
    if quiet:
        loglevel = "WARNING"
    log.init_log(logger, loglevel=loglevel, no_ansi=no_ansi)

    console = rich.console.Console(stderr=True, highlight=False, force_terminal=log.force_term_colors())
    console.print(
        "\n  [blue]/[/][green]/[/][red]/[/] [bold][link=https://multiqc.info]MultiQC[/link][/] :mag: [dim]| v{}\n".format(
            config.version
        )
    )
    logger.debug("This is MultiQC v{}".format(config.version))

    # Load config files
    plugin_hooks.mqc_trigger("before_config")
    config.mqc_load_userconfig(config_file)
    plugin_hooks.mqc_trigger("config_loaded")

    # Command-line config YAML
    if len(cl_config) > 0:
        config.mqc_cl_config(cl_config)

    # Log the command used to launch MultiQC
    report.multiqc_command = " ".join(sys.argv)
    logger.debug("Command used: {}".format(report.multiqc_command))

    # Check that we're running the latest version of MultiQC
    if config.no_version_check is not True:
        try:
            response = urlopen("http://multiqc.info/version.php?v={}".format(config.short_version), timeout=5)
            remote_version = response.read().decode("utf-8").strip()
            if version.StrictVersion(re.sub("[^0-9\.]", "", remote_version)) > version.StrictVersion(
                re.sub("[^0-9\.]", "", config.short_version)
            ):
                logger.warning("MultiQC Version {} now available!".format(remote_version))
            else:
                logger.debug("Latest MultiQC version is {}".format(remote_version))
        except Exception as e:
            logger.debug("Could not connect to multiqc.info for version check: {}".format(e))

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
    config.analysis_dir = analysis_dir
    if outdir is not None:
        config.output_dir = outdir
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
    if plots_flat:
        config.plots_force_flat = True
    if plots_interactive:
        config.plots_force_interactive = True
    if lint:
        config.lint = True
        lint_helpers.run_tests()
    if make_pdf:
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
    if module_tag is not None:
        config.module_tag = module_tag
    if len(module) > 0:
        config.run_modules = module
    if len(exclude) > 0:
        config.exclude_modules = exclude
    if profile_runtime:
        config.profile_runtime = True
    config.kwargs = kwargs  # Plugin command line options

    # Clean up analysis_dir if a string (interactive environment only)
    if isinstance(config.analysis_dir, str):
        config.analysis_dir = [config.analysis_dir]

    plugin_hooks.mqc_trigger("execution_start")

    logger.debug("Command     : {}".format(" ".join(sys.argv)))
    logger.debug("Working dir : {}".format(os.getcwd()))
    if make_pdf:
        logger.info("--pdf specified. Using non-interactive HTML template.")
    logger.debug("Template    : {}".format(config.template))
    if lint:
        logger.info("--lint specified. Being strict with validation.")

    # Throw a warning if we are running on Python 2
    if sys.version_info[0] < 3:
        logger.error(
            "You are running MultiQC with Python {}.{}.{}".format(
                sys.version_info[0], sys.version_info[1], sys.version_info[2]
            )
        )
        logger.critical("Please upgrade Python! MultiQC does not support Python < 3.6, things will break.")
    else:
        logger.debug("Running Python {}".format(sys.version.replace("\n", " ")))

    # Add files if --file-list option is given
    if file_list:
        if len(analysis_dir) > 1:
            raise ValueError("If --file-list is giving, analysis_dir should have only one plain text file.")
        config.analysis_dir = []
        with (open(analysis_dir[0])) as in_handle:
            for line in in_handle:
                if os.path.exists(line.strip()):
                    path = os.path.abspath(line.strip())
                    config.analysis_dir.append(path)
        if len(config.analysis_dir) == 0:
            logger.error("No files or directories were added from {} using --file-list option.".format(analysis_dir[0]))
            logger.error("Please, check that {} contains correct paths.".format(analysis_dir[0]))
            raise ValueError("Any files or directories to be searched.")

    if len(ignore) > 0:
        logger.debug("Ignoring files, directories and paths that match: {}".format(", ".join(ignore)))
        config.fn_ignore_files.extend(ignore)
        config.fn_ignore_dirs.extend(ignore)
        config.fn_ignore_paths.extend(ignore)
    if len(ignore_samples) > 0:
        logger.debug("Ignoring sample names that match: {}".format(", ".join(ignore_samples)))
        config.sample_names_ignore.extend(ignore_samples)
    if filename == "stdout":
        config.output_fn = sys.stdout
        logger.info("Printing report to stdout")
    else:
        if title is not None and filename is None:
            filename = re.sub("[^\w\.-]", "", re.sub("[-\s]+", "-", title)).strip()
            filename += "_multiqc_report"
        if filename is not None:
            if filename.endswith(".html"):
                filename = filename[:-5]
            config.output_fn_name = filename
            config.data_dir_name = "{}_data".format(filename)
            config.plots_dir_name = "{}_plots".format(filename)
        if not config.output_fn_name.endswith(".html"):
            config.output_fn_name = "{}.html".format(config.output_fn_name)

    # Print some status updates
    if config.title is not None:
        logger.info("Report title: {}".format(config.title))
    if dirs:
        logger.info("Prepending directory to sample names")

    # Prep module configs
    config.top_modules = [m if type(m) is dict else {m: {}} for m in config.top_modules]
    config.module_order = [m if type(m) is dict else {m: {}} for m in config.module_order]
    mod_keys = [list(m.keys())[0] for m in config.module_order]

    # Lint the module configs
    if config.lint:
        for m in config.avail_modules.keys():
            if m not in mod_keys:
                errmsg = "LINT: Module '{}' not found in config.module_order".format(m)
                logger.error(errmsg)
                report.lint_errors.append(errmsg)
            else:
                for mo in config.module_order:
                    if m != "custom_content" and m in mo.keys() and "module_tag" not in mo[m]:
                        errmsg = "LINT: Module '{}' in config.module_order did not have 'module_tag' config".format(m)
                        logger.error(errmsg)
                        report.lint_errors.append(errmsg)

    # Get the avaiable tags to decide which modules to run.
    modules_from_tags = set()
    if config.module_tag is not None:
        tags = config.module_tag
        for m in config.module_order:
            module_name = list(m.keys())[0]  # only one name in each dict
            for tag in tags:
                for t in m[module_name].get("module_tag", []):
                    if tag.lower() == t.lower():
                        modules_from_tags.add(module_name)

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

    if len(getattr(config, "run_modules", {})) > 0:
        run_modules = [m for m in run_modules if list(m.keys())[0] in config.run_modules]
        logger.info("Only using modules: {}".format(", ".join(config.run_modules)))
    elif modules_from_tags:
        run_modules = [m for m in run_modules if list(m.keys())[0] in modules_from_tags]
        logger.info("Only using modules with '{}' tag".format(", ".join(module_tag)))
    if len(getattr(config, "exclude_modules", {})) > 0:
        logger.info("Excluding modules '{}'".format("', '".join(config.exclude_modules)))
        if "general_stats" in config.exclude_modules:
            config.skip_generalstats = True
            config.exclude_modules = tuple(x for x in config.exclude_modules if x != "general_stats")
        run_modules = [m for m in run_modules if list(m.keys())[0] not in config.exclude_modules]
    if len(run_modules) == 0:
        logger.critical("No analysis modules specified!")
        sys.exit(1)
    run_module_names = [list(m.keys())[0] for m in run_modules]
    logger.debug("Analysing modules: {}".format(", ".join(run_module_names)))

    # Create the temporary working directories
    tmp_dir = tempfile.mkdtemp()
    logger.debug("Using temporary directory for creating report: {}".format(tmp_dir))
    config.data_tmp_dir = os.path.join(tmp_dir, "multiqc_data")
    if filename != "stdout" and config.make_data_dir == True:
        config.data_dir = config.data_tmp_dir
        os.makedirs(config.data_dir)
    else:
        config.data_dir = None
    config.plots_tmp_dir = os.path.join(tmp_dir, "multiqc_plots")
    if filename != "stdout" and config.export_plots == True:
        config.plots_dir = config.plots_tmp_dir
        os.makedirs(config.plots_dir)
    else:
        config.plots_dir = None

    # Load the template
    template_mod = config.avail_templates[config.template].load()

    # Add an output subdirectory if specified by template
    try:
        config.output_dir = os.path.join(config.output_dir, template_mod.output_subdir)
    except AttributeError:
        pass  # No subdirectory variable given

    # Add custom content section names
    try:
        if "custom_content" in run_module_names:
            run_module_names.extend(config.custom_data.keys())
    except AttributeError:
        pass  # custom_data not in config

    # Get the list of files to search
    for d in config.analysis_dir:
        logger.info("Search path : {}".format(os.path.abspath(d)))
    report.get_filelist(run_module_names)

    # Run the modules!
    plugin_hooks.mqc_trigger("before_modules")
    report.modules_output = list()
    sys_exit_code = 0
    total_mods_starttime = time.time()
    for mod_idx, mod_dict in enumerate(run_modules):
        mod_starttime = time.time()
        try:
            this_module = list(mod_dict.keys())[0]
            mod_cust_config = list(mod_dict.values())[0]
            if mod_cust_config is None:
                mod_cust_config = {}
            mod = config.avail_modules[this_module].load()
            mod.mod_cust_config = mod_cust_config  # feels bad doing this, but seems to work
            output = mod()
            if type(output) != list:
                output = [output]
            for m in output:
                report.modules_output.append(m)

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

        except UserWarning:
            logger.debug("No samples found: {}".format(list(mod_dict.keys())[0]))
        except KeyboardInterrupt:
            shutil.rmtree(tmp_dir)
            logger.critical(
                "User Cancelled Execution!\n{eq}\n{tb}{eq}\n".format(eq=("=" * 60), tb=traceback.format_exc())
                + "User Cancelled Execution!\nExiting MultiQC..."
            )
            sys.exit(1)
        except:
            # Flag the error, but carry on
            class CustomTraceback:
                def __rich_console__(self, console: rich.console.Console, options: rich.console.ConsoleOptions):
                    sys_tb = sys.exc_info()
                    issue_url = "https://github.com/ewels/MultiQC/issues/new?template=bug_report.md&title={}%20module%20-%20{}".format(
                        this_module, sys_tb[0].__name__
                    )
                    yield (
                        "Please copy this log and report it at [bright_blue][link={}]https://github.com/ewels/MultiQC/issues[/link][/] \n"
                        "[bold underline]Please attach a file that triggers the error.[/] The last file found was: [green]{}[/]\n".format(
                            issue_url, report.last_found_file
                        )
                    )
                    yield Syntax(traceback.format_exc(), "python")

                def __rich_measure__(self, console: rich.console.Console, options: rich.console.ConsoleOptions):
                    tb_width = max([len(l) for l in traceback.format_exc().split("\n")])
                    try:
                        log_width = 71 + len(report.last_found_file)
                    except TypeError:
                        log_width = 71
                    panel_width = max(tb_width, log_width)
                    return rich.console.Measurement(panel_width, panel_width)

            console = rich.console.Console(stderr=True, force_terminal=log.force_term_colors())
            console.print(
                rich.panel.Panel(
                    CustomTraceback(),
                    title="Oops! The '[underline]{}[/]' MultiQC module broke...".format(this_module),
                    expand=False,
                    border_style="red",
                    style="on #272822",
                )
            )
            # Still log.debug this so that it ends up in the log file - above is just stderr for now
            logger.debug(
                "Oops! The '{}' MultiQC module broke...\n".format(this_module)
                + ("=" * 80)
                + "\n"
                + traceback.format_exc()
                + ("=" * 80)
            )
            # Exit code 1 for CI failures etc
            sys_exit_code = 1

        report.runtimes["mods"][run_module_names[mod_idx]] = time.time() - mod_starttime
    report.runtimes["total_mods"] = time.time() - total_mods_starttime

    # Special-case module if we want to profile the MultiQC running time
    if config.profile_runtime:
        from multiqc.utils import profile_runtime

        report.modules_output.append(profile_runtime.MultiqcModule())

    # Did we find anything?
    if len(report.modules_output) == 0:
        logger.warning("No analysis results found. Cleaning up..")
        shutil.rmtree(tmp_dir)
        logger.info("MultiQC complete")
        # Exit with an error code if a module broke
        sys.exit(sys_exit_code)

    # Sort the report module output if we have a config
    if len(getattr(config, "report_section_order", {})) > 0:
        section_id_order = {}
        idx = 10
        for mod in reversed(report.modules_output):
            section_id_order[mod.anchor] = idx
            idx += 10
        for anchor, ss in config.report_section_order.items():
            if anchor not in section_id_order.keys():
                logger.debug("Reordering sections: anchor '{}' not found.".format(anchor))
                continue
            if ss.get("order") is not None:
                section_id_order[anchor] = ss["order"]
            if ss.get("after") in section_id_order.keys():
                section_id_order[anchor] = section_id_order[ss["after"]] + 1
            if ss.get("before") in section_id_order.keys():
                section_id_order[anchor] = section_id_order[ss["before"]] - 1
        sorted_ids = sorted(section_id_order, key=section_id_order.get)
        report.modules_output = [mod for i in reversed(sorted_ids) for mod in report.modules_output if mod.anchor == i]

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
                    logger.debug("Reordering sections: anchor '{}' not found for module '{}'.".format(anchor, mod.name))
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
                "mqc-generalstats-{}-{}".format(ns_html, h[k]["rid"])
            )
    # Generate the General Statistics HTML & write to file
    if len(report.general_stats_data) > 0 and not config.skip_generalstats:
        pconfig = {
            "id": "general_stats_table",
            "table_title": "General Statistics",
            "save_file": True,
            "raw_data_fn": "multiqc_general_stats",
        }
        report.general_stats_html = table.plot(report.general_stats_data, report.general_stats_headers, pconfig)
    else:
        config.skip_generalstats = True

    # Write the report sources to disk
    if config.data_dir is not None:
        report.data_sources_tofile()
    # Compress the report plot JSON data
    runtime_compression_start = time.time()
    logger.info("Compressing plot data")
    report.plot_compressed_json = report.compress_json(report.plot_data)
    report.runtimes["total_compression"] = time.time() - runtime_compression_start

    plugin_hooks.mqc_trigger("before_report_generation")

    # Data Export / MegaQC integration - save report data to file or send report data to an API endpoint
    if (config.data_dump_file or config.megaqc_url) and config.megaqc_upload:
        multiqc_json_dump = megaqc.multiqc_dump_json(report)
        if config.data_dump_file:
            util_functions.write_data_file(multiqc_json_dump, "multiqc_data", False, "json")
        if config.megaqc_url:
            megaqc.multiqc_api_post(multiqc_json_dump)

    # Make the final report path & data directories
    if filename != "stdout":
        config.output_fn = os.path.join(config.output_dir, config.output_fn_name)
        config.data_dir = os.path.join(config.output_dir, config.data_dir_name)
        config.plots_dir = os.path.join(config.output_dir, config.plots_dir_name)
        # Check for existing reports and remove if -f was specified
        if (
            os.path.exists(config.output_fn)
            or (config.make_data_dir and os.path.exists(config.data_dir))
            or (config.export_plots and os.path.exists(config.plots_dir))
        ):
            if config.force:
                if os.path.exists(config.output_fn):
                    logger.warning("Deleting    : {}   (-f was specified)".format(os.path.relpath(config.output_fn)))
                    os.remove(config.output_fn)
                if config.make_data_dir and os.path.exists(config.data_dir):
                    logger.warning("Deleting    : {}   (-f was specified)".format(os.path.relpath(config.data_dir)))
                    shutil.rmtree(config.data_dir)
                if config.export_plots and os.path.exists(config.plots_dir):
                    logger.warning("Deleting    : {}   (-f was specified)".format(os.path.relpath(config.plots_dir)))
                    shutil.rmtree(config.plots_dir)
            else:
                # Set up the base names of the report and the data dir
                report_num = 1
                report_base, report_ext = os.path.splitext(config.output_fn_name)
                dir_base = os.path.basename(config.data_dir)
                plots_base = os.path.basename(config.plots_dir)

                # Iterate through appended numbers until we find one that's free
                while (
                    os.path.exists(config.output_fn)
                    or (config.make_data_dir and os.path.exists(config.data_dir))
                    or (config.export_plots and os.path.exists(config.plots_dir))
                ):
                    config.output_fn = os.path.join(
                        config.output_dir, "{}_{}{}".format(report_base, report_num, report_ext)
                    )
                    config.data_dir = os.path.join(config.output_dir, "{}_{}".format(dir_base, report_num))
                    config.plots_dir = os.path.join(config.output_dir, "{}_{}".format(plots_base, report_num))
                    report_num += 1

                config.output_fn_name = os.path.basename(config.output_fn)
                config.data_dir_name = os.path.basename(config.data_dir)
                config.plots_dir_name = os.path.basename(config.plots_dir)
                logger.warning("Previous MultiQC output found! Adjusting filenames..")
                logger.warning("Use -f or --force to overwrite existing reports instead")

        # Make directories for report if needed
        if not os.path.exists(os.path.dirname(config.output_fn)):
            os.makedirs(os.path.dirname(config.output_fn))
        logger.info("Report      : {}".format(os.path.relpath(config.output_fn)))

        if config.make_data_dir == False:
            logger.info("Data        : None")
        else:
            # Make directories for data_dir
            logger.info("Data        : {}".format(os.path.relpath(config.data_dir)))
            # Modules have run, so data directory should be complete by now. Move its contents.
            logger.debug("Moving data file from '{}' to '{}'".format(config.data_tmp_dir, config.data_dir))
            # Disable preserving of times and mode on purpose to avoid problems with mounted CIFS shares (see #625)
            copy_tree(config.data_tmp_dir, config.data_dir, preserve_times=0, preserve_mode=0)
            shutil.rmtree(config.data_tmp_dir)

        # Copy across the static plot images if requested
        if config.export_plots:
            config.plots_dir = os.path.join(config.output_dir, config.plots_dir_name)
            if os.path.exists(config.plots_dir):
                if config.force:
                    logger.warning("Deleting    : {}   (-f was specified)".format(os.path.relpath(config.plots_dir)))
                    shutil.rmtree(config.plots_dir)
                else:
                    logger.error("Output directory {} already exists.".format(config.plots_dir))
                    logger.info("Use -f or --force to overwrite existing reports")
                    shutil.rmtree(tmp_dir)
                    sys.exit(1)
            logger.info("Plots       : {}".format(os.path.relpath(config.plots_dir)))

            # Modules have run, so plots directory should be complete by now. Move its contents.
            logger.debug("Moving plots directory from '{}' to '{}'".format(config.plots_tmp_dir, config.plots_dir))
            # Disable preserving of times and mode on purpose to avoid problems with mounted CIFS shares (see #625)
            copy_tree(config.plots_tmp_dir, config.plots_dir, preserve_times=0, preserve_mode=0)
            shutil.rmtree(config.plots_tmp_dir)

    plugin_hooks.mqc_trigger("before_template")

    # Load in parent template files first if a child theme
    try:
        parent_template = config.avail_templates[template_mod.template_parent].load()
        copy_tree(parent_template.template_dir, tmp_dir)
    except AttributeError:
        pass  # Not a child theme

    # Copy the template files to the tmp directory (distutils overwrites parent theme files)
    copy_tree(template_mod.template_dir, tmp_dir)

    # Function to include file contents in Jinja template
    def include_file(name, fdir=tmp_dir, b64=False):
        try:
            if fdir is None:
                fdir = ""
            if b64:
                with io.open(os.path.join(fdir, name), "rb") as f:
                    return base64.b64encode(f.read()).decode("utf-8")
            else:
                with io.open(os.path.join(fdir, name), "r", encoding="utf-8") as f:
                    return f.read()
        except (OSError, IOError) as e:
            logger.error("Could not include file '{}': {}".format(name, e))

    # Load the report template
    try:
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(tmp_dir))
        env.globals["include_file"] = include_file
        j_template = env.get_template(template_mod.base_fn)
    except:
        raise IOError("Could not load {} template file '{}'".format(config.template, template_mod.base_fn))

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
            raise IOError("Could not print report to '{}' - {}".format(config.output_fn, IOError(e)))

        # Copy over files if requested by the theme
        try:
            for f in template_mod.copy_files:
                fn = os.path.join(tmp_dir, f)
                dest_dir = os.path.join(os.path.dirname(config.output_fn), f)
                copy_tree(fn, dest_dir)
        except AttributeError:
            pass  # No files to copy

    # Clean up temporary directory
    shutil.rmtree(tmp_dir)

    # Zip the data directory if requested
    if config.zip_data_dir and config.data_dir is not None:
        shutil.make_archive(config.data_dir, "zip", config.data_dir)
        shutil.rmtree(config.data_dir)

    # Try to create a PDF if requested
    if make_pdf:
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
                pandoc_call.append("--template={}".format(config.pandoc_template))
            logger.debug(
                "Attempting Pandoc conversion to PDF with following command:\n{}".format(" ".join(pandoc_call))
            )
            pdf_exit_code = subprocess.call(pandoc_call)
            if pdf_exit_code != 0:
                logger.error("Error creating PDF! Pandoc returned a non-zero exit code.")
            else:
                logger.info("PDF Report  : {}".format(pdf_fn_name))
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                logger.error("Error creating PDF - pandoc not found. Is it installed? http://pandoc.org/")
            else:
                logger.error(
                    "Error creating PDF! Something went wrong when creating the PDF\n"
                    + ("=" * 60)
                    + "\n{}\n".format(traceback.format_exc())
                    + ("=" * 60)
                )

    plugin_hooks.mqc_trigger("execution_finish")

    logger.info("MultiQC complete")
    report.runtimes["total"] = time.time() - start_execution_time
    if config.profile_runtime:
        logger.info("Run took {:.2f} seconds".format(report.runtimes["total"]))
        logger.info(" - {:.2f}s: Searching files".format(report.runtimes["total_sp"]))
        logger.info(" - {:.2f}s: Running modules".format(report.runtimes["total_mods"]))
        logger.info(" - {:.2f}s: Compressing report data".format(report.runtimes["total_compression"]))
        logger.info("For more information, see the 'Run Time' section in {}".format(os.path.relpath(config.output_fn)))

    if report.num_mpl_plots > 0 and not config.plots_force_flat:
        logger.warning(
            "{} flat-image plot{} used in the report due to large sample numbers".format(
                report.num_mpl_plots, "s" if report.num_mpl_plots > 1 else ""
            )
        )
        console.print(
            "[blue]|           multiqc[/] | "
            "To force interactive plots, use the [yellow]'--interactive'[/] flag. "
            "See the [link=https://multiqc.info/docs/#flat--interactive-plots]documentation[/link]."
        )

    if lint and len(report.lint_errors) > 0:
        logger.error("Found {} linting errors!\n{}".format(len(report.lint_errors), "\n".join(report.lint_errors)))
        sys_exit_code = 1

    # Move the log file into the data directory
    log.move_tmp_log(logger)

    # Return the running information from the run:
    #
    # * report instance
    # * config instance
    # * appropriate error code (eg. 1 if a module broke, 0 on success)
    #
    return {"report": report, "config": config, "sys_exit_code": sys_exit_code}
