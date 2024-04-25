import logging
import os
import platform
import re
import sys

import requests
from packaging import version

from multiqc.core.exceptions import RunError
from multiqc.utils import report, config, plugin_hooks, strict_helpers
from multiqc.utils.util_functions import strtobool
from multiqc.utils import log

logger = logging.getLogger(__name__)


def init_config(
    analysis_dir=None,
    file_list=None,
    prepend_dirs=None,
    dirs_depth=None,
    fn_clean_sample_names=None,
    title=None,
    report_comment=None,
    template=None,
    run_modules=(),
    exclude_modules=(),
    require_logs=None,
    output_dir=None,
    use_filename_as_sample_name=False,
    replace_names=None,
    sample_names=None,
    sample_filters=None,
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
    development=None,
    make_pdf=None,
    no_megaqc_upload=None,
    config_file=(),
    cl_config=(),
    quiet=None,
    verbose=None,
    no_ansi=None,
    profile_runtime=None,
    custom_css_files=(),
    module_order=(),
    **kwargs,
):
    """
    Initialize log and config.
    Update config from command line parameters. Will only update from the non-None values.
    """
    config.load_from_defaults()

    if quiet is not None:
        config.quiet = quiet
    if no_ansi is not None:
        config.no_ansi = no_ansi
    if verbose is not None:
        config.verbose = verbose > 0
    log.init_log()
    logger.debug(f"This is MultiQC v{config.version}")

    plugin_hooks.mqc_trigger("before_config")
    config.load_userconfig(config_file)
    plugin_hooks.mqc_trigger("config_loaded")

    # Command-line config YAML
    if len(cl_config) > 0:
        config._cl_config(cl_config)

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
        logger.info(f"Report title: {config.title}")
    if report_comment is not None:
        config.report_comment = report_comment
    if prepend_dirs is not None:
        config.prepend_dirs = prepend_dirs
        logger.info("Prepending directory to sample names")
    if dirs_depth is not None:
        config.prepend_dirs = True
        config.prepend_dirs_depth = dirs_depth
    if output_dir is not None:
        config.output_dir = os.path.realpath(output_dir)
    if use_filename_as_sample_name is not None:
        config.use_filename_as_sample_name = use_filename_as_sample_name
        logger.info("Using log filenames for sample names")
    if make_data_dir is not None:
        config.make_data_dir = make_data_dir
    if force is not None:
        config.force = force
    if ignore_symlinks is not None:
        config.ignore_symlinks = ignore_symlinks
    if zip_data_dir is not None:
        config.zip_data_dir = zip_data_dir
    if data_format is not None:
        config.data_format = data_format
    if export_plots is not None:
        config.export_plots = export_plots
    if make_report is not None:
        config.make_report = make_report
    if plots_force_flat is not None:
        config.plots_force_flat = plots_force_flat
    if plots_force_interactive is not None:
        config.plots_force_interactive = plots_force_interactive
    if strict is not None:
        config.strict = strict
        config.lint = strict  # Deprecated since v1.17
        if strict:
            strict_helpers.run_tests()
    if development is not None:
        config.development = development
        if development:
            if "png" not in config.export_plot_formats:
                config.export_plot_formats.append("png")
    if make_pdf:
        config.template = "simple"
    if filename:
        config.filename = filename
    if no_megaqc_upload is not None:
        config.megaqc_upload = not no_megaqc_upload
    if fn_clean_sample_names is not None:
        config.fn_clean_sample_names = fn_clean_sample_names
        logger.info("Not cleaning sample names")
    if replace_names:
        config.load_replace_names(replace_names)
    if sample_names:
        config.load_sample_names(sample_names)
    config.load_show_hide(sample_filters)
    if len(run_modules) > 0:
        config.run_modules = run_modules
    if len(exclude_modules) > 0:
        config.exclude_modules = exclude_modules
    if require_logs is not None:
        config.require_logs = require_logs
    if profile_runtime is not None:
        config.profile_runtime = profile_runtime
    if custom_css_files:
        config.custom_css_files.extend(custom_css_files)
    if module_order:
        config.module_order = module_order

    _set_analysis_file_config(
        analysis_dir=analysis_dir,
        file_list=file_list,
        ignore=(),
        ignore_samples=(),
    )

    _set_output_paths()

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


def _set_analysis_file_config(
    analysis_dir=None,
    file_list=None,
    ignore=(),
    ignore_samples=(),
):
    # Clean up analysis_dir if a string (interactive environment only)
    if analysis_dir is not None:
        if isinstance(analysis_dir, str):
            analysis_dir = [analysis_dir]
        config.analysis_dir = analysis_dir
    if file_list is not None:
        config.file_list = file_list

    # Add files if --file-list option is given
    if config.file_list:
        if len(config.analysis_dir) > 1:
            raise RunError("If --file-list is given, analysis_dir should have only one plain text file.")
        file_list_path = config.analysis_dir[0]
        config.analysis_dir = []
        with open(file_list_path) as in_handle:
            for line in in_handle:
                if os.path.exists(line.strip()):
                    path = os.path.abspath(line.strip())
                    config.analysis_dir.append(path)
        if len(config.analysis_dir) == 0:
            raise RunError(
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

    # Prep module configs
    config.top_modules = [m if isinstance(m, dict) else {m: {}} for m in config.top_modules]
    config.module_order = [m if isinstance(m, dict) else {m: {}} for m in config.module_order]

    # Lint the module config
    mod_keys = [list(m.keys())[0] for m in config.module_order]
    if config.strict:
        for m in config.avail_modules.keys():
            if m not in mod_keys:
                errmsg = f"LINT: Module '{m}' not found in config.module_order"
                logger.error(errmsg)
                report.lint_errors.append(errmsg)


def _set_output_paths():
    """
    Create directories where outputs will be placed
    """

    # Add an output subdirectory if specified by template
    template_mod = config.avail_templates[config.template].load()
    try:
        config.output_dir = os.path.join(config.output_dir, template_mod.output_subdir)
    except AttributeError:
        pass  # No subdirectory variable given

    filename = config.filename

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
        else:
            config.output_fn = None
        config.data_dir = os.path.join(config.output_dir, config.data_dir_name)
        config.plots_dir = os.path.join(config.output_dir, config.plots_dir_name)
