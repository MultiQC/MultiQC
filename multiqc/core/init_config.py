import logging
import os
import sys
from pathlib import Path

from multiqc import report, config
from multiqc.core.exceptions import RunError
from multiqc.core import init_log, plugin_hooks, strict_helpers

logger = logging.getLogger(__name__)


def update_config(
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
    ignore=(),
    ignore_samples=(),
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
    config_files=(),
    cl_config=(),
    quiet=None,
    verbose=None,
    no_ansi=None,
    profile_runtime=None,
    custom_css_files=(),
    module_order=(),
    extra_fn_clean_exts=(),
    extra_fn_clean_trim=(),
    **kwargs,
):
    """
    Update config and logger from parameters.

    Will reload config from defaults and from the previously added user config files,
    and then update from the non-None arguments above.
    """

    # Reload from defaults
    config.load_defaults()

    # Set up logger
    if quiet is not None:
        config.quiet = quiet
    if no_ansi is not None:
        config.no_ansi = no_ansi
    if verbose is not None:
        config.verbose = verbose > 0
    init_log.init_log()
    logger.debug(f"This is MultiQC v{config.version}")
    logger.debug("Running Python " + sys.version.replace("\n", " "))

    plugin_hooks.mqc_trigger("before_config")

    # Set up user configs (they are kept for the whole interactive session)
    for path in config_files:
        if path not in config.session_user_config_files:
            config.session_user_config_files.append(Path(path).absolute())

    # Command-line config YAML
    if len(cl_config) > 0:
        config.load_cl_config(cl_config)

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
        if use_filename_as_sample_name:
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
    if extra_fn_clean_exts:
        config.fn_clean_exts = list(extra_fn_clean_exts) + config.fn_clean_exts
    if extra_fn_clean_trim:
        config.fn_clean_trim = list(extra_fn_clean_trim) + config.fn_clean_trim

    # Clean up analysis_dir if a string (interactive environment only)
    if analysis_dir is not None:
        if isinstance(analysis_dir, str):
            analysis_dir = [analysis_dir]
        config.analysis_dir = analysis_dir
    if file_list is not None:
        if len(config.analysis_dir) > 1:
            raise RunError("If --file-list is given, analysis_dir should have only one plain text file.")
        config.file_list = file_list

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

    config.kwargs = kwargs  # Plugin command line options

    plugin_hooks.mqc_trigger("config_loaded")
    plugin_hooks.mqc_trigger("execution_start")
