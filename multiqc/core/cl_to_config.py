import logging
import os
import platform
import re
import sys

import requests
from packaging import version

from multiqc.utils import report, config, plugin_hooks, strict_helpers
from multiqc.utils.util_functions import strtobool

logger = logging.getLogger(__name__)


def _cl_to_config(
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
    profile_runtime=False,
    no_ansi=False,
    custom_css_files=(),
    module_order=(),
    **kwargs,
):
    """
    Populate config from command line parameters.
    """
    plugin_hooks.mqc_trigger("before_config")
    config.mqc_load_userconfig(config_file)
    plugin_hooks.mqc_trigger("config_loaded")

    # Command-line config YAML
    if len(cl_config) > 0:
        config.mqc_cl_config(cl_config)

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
        if "png" not in config.export_plot_formats:
            config.export_plot_formats.append("png")
    if make_pdf:
        config.template = "simple"
    if filename:
        config.filename = filename
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
    if module_order:
        config.module_order = module_order

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
