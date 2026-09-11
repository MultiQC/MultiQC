import logging
import os
import sys
from pathlib import Path
from typing import Literal, cast

from pydantic import BaseModel

from multiqc import config, report
from multiqc.core import log_and_rich, plugin_hooks
from multiqc.core.exceptions import RunError
from multiqc.utils.config_schema import AiProviderLiteral

logger = logging.getLogger(__name__)


class ClConfig(BaseModel):
    """
    Holds config updates from the command line or interactive functions.
    """

    file_list: bool | None = None
    prepend_dirs: bool | None = None
    dirs_depth: int | None = None
    fn_clean_sample_names: bool | None = None
    title: str | None = None
    report_comment: str | None = None
    template: str | None = None
    require_logs: bool | None = None
    output_dir: str | Path | None = None
    use_filename_as_sample_name: bool | list[str] | None = None
    replace_names: str | None = None
    sample_names: str | None = None
    sample_filters: str | None = None
    filename: str | None = None
    make_data_dir: bool | None = None
    data_format: str | None = None
    zip_data_dir: bool | None = None
    force: bool | None = None
    ignore_symlinks: bool | None = None
    make_report: bool | None = None
    export_plots: bool | None = None
    plots_force_flat: bool | None = None
    plots_force_interactive: bool | None = None
    strict: bool | None = None
    development: bool | None = None
    make_pdf: bool | None = None
    no_megaqc_upload: bool | None = None
    quiet: bool | None = None
    verbose: bool | None = None
    no_ansi: bool | None = None
    profile_runtime: bool | None = None
    profile_memory: bool | None = None
    no_version_check: bool | None = None
    ignore: list[str] = []
    ignore_samples: list[str] = []
    only_samples: list[str] = []
    run_modules: list[str] = []
    exclude_modules: list[str] = []
    config_files: list[str | Path] = []
    cl_config: list[str] = []
    custom_css_files: list[str] = []
    module_order: list[str | dict] = []
    extra_fn_clean_exts: list = []
    extra_fn_clean_trim: list = []
    preserve_module_raw_data: bool | None = None
    data_dump_file_write_raw: bool | None = None
    ai_summary: bool | None = None
    ai_summary_full: bool | None = None
    ai_provider: AiProviderLiteral | None = None
    ai_model: str | None = None
    ai_custom_endpoint: str | None = None
    ai_custom_context_window: int | None = None
    ai_prompt_short: str | None = None
    ai_prompt_full: str | None = None
    no_ai: bool | None = None
    unknown_options: dict | None = None
    check_config: bool | None = None


def update_config(*analysis_dir, cfg: ClConfig | None = None, log_to_file=False, print_intro_fn=None):
    """
    Update config and re-initialize logger.

    First will reload config from defaults and from the previously added user config files.
    Then will update from cfg from the non-None arguments.
    """

    # Reload from defaults
    config.load_defaults()

    if cfg is None and analysis_dir and isinstance(analysis_dir[0], ClConfig):
        cfg = analysis_dir[0]
    cfg = cfg or ClConfig()

    # Reset logger
    if cfg.quiet is not None:
        config.quiet = cfg.quiet
    if cfg.no_ansi is not None:
        config.no_ansi = cfg.no_ansi
    if cfg.verbose is not None:
        config.verbose = cfg.verbose > 0

    log_and_rich.init_log(log_to_file=log_to_file)
    if print_intro_fn is not None:
        print_intro_fn()

    logger.debug(f"This is MultiQC v{config.version}")
    logger.debug("Running Python " + sys.version.replace("\n", " "))

    plugin_hooks.mqc_trigger("before_config")

    config.loaded_user_files = set()

    # Re-finding implicit configs
    config.find_user_files()

    # Re-loading explicit user configs
    path: Path | str
    for path in config.explicit_user_config_files:
        config.load_config_file(path)

    # Set up session config files passed with -c or cfg=
    for path in cfg.config_files:
        config.load_config_file(str(path), is_explicit_config=False)

    # Command-line config YAML
    if len(cfg.cl_config) > 0:
        config.load_cl_config(cfg.cl_config)

    # Set up key variables (overwrite config vars from command line)
    if cfg.template is not None:
        # `cfg.template` is a `str` from click; click.Choice already
        # validated it against the available template names.
        config.template = cast(
            Literal["default", "original", "simple", "sections", "gathered", "geo", "disco"],
            cfg.template,
        )
    if cfg.title is not None:
        config.title = cfg.title
        logger.info(f"Report title: {config.title}")
    if cfg.report_comment is not None:
        config.report_comment = cfg.report_comment
    if cfg.prepend_dirs is not None:
        config.prepend_dirs = cfg.prepend_dirs
        if cfg.prepend_dirs:
            logger.info("Prepending directory to sample names")
    if cfg.dirs_depth is not None:
        config.prepend_dirs = True
        config.prepend_dirs_depth = cfg.dirs_depth
    if cfg.output_dir is not None:
        config.output_dir = os.path.realpath(cfg.output_dir)
    if cfg.use_filename_as_sample_name is not None:
        config.use_filename_as_sample_name = cfg.use_filename_as_sample_name
        if cfg.use_filename_as_sample_name:
            logger.info("Using log filenames for sample names")
    if cfg.make_data_dir is not None:
        config.make_data_dir = cfg.make_data_dir
    if cfg.force is not None:
        config.force = cfg.force
    if cfg.ignore_symlinks is not None:
        config.ignore_symlinks = cfg.ignore_symlinks
    if cfg.zip_data_dir is not None:
        config.zip_data_dir = cfg.zip_data_dir
    if cfg.data_format is not None:
        config.data_format = cast(Literal["tsv", "csv", "json", "yaml"], cfg.data_format)
    if cfg.export_plots is not None:
        config.export_plots = cfg.export_plots
    if cfg.make_report is not None:
        config.make_report = cfg.make_report
    if cfg.plots_force_flat is not None:
        config.plots_force_flat = cfg.plots_force_flat
    if cfg.plots_force_interactive is not None:
        config.plots_force_interactive = cfg.plots_force_interactive
    if cfg.strict is not None:
        config.strict = cfg.strict
        config.lint = cfg.strict  # Deprecated since v1.17
    if cfg.development is not None:
        config.development = cfg.development
    if cfg.make_pdf:
        config.make_pdf = cfg.make_pdf
        config.template = "simple"
    if config.template == "simple":
        config.plots_force_flat = True
        config.simple_output = True
    if cfg.filename:
        config.filename = cfg.filename
    if cfg.no_megaqc_upload is not None:
        config.megaqc_upload = not cfg.no_megaqc_upload
    if cfg.fn_clean_sample_names is not None:
        config.fn_clean_sample_names = cfg.fn_clean_sample_names
        if not cfg.fn_clean_sample_names:
            logger.info("Not cleaning sample names")
    if cfg.replace_names:
        config.load_replace_names(Path(cfg.replace_names))
    if cfg.sample_names:
        config.load_sample_names(Path(cfg.sample_names))
    config.load_show_hide(show_hide_file=Path(cfg.sample_filters) if cfg.sample_filters else None)
    if len(cfg.run_modules) > 0:
        config.run_modules = cfg.run_modules
    if len(cfg.exclude_modules) > 0:
        config.exclude_modules = cfg.exclude_modules
    if cfg.require_logs is not None:
        config.require_logs = cfg.require_logs
    if cfg.profile_runtime is not None:
        config.profile_runtime = cfg.profile_runtime
    if cfg.profile_memory is not None:
        config.profile_runtime = config.profile_memory = cfg.profile_memory
    if cfg.no_version_check is not None:
        config.no_version_check = cfg.no_version_check
    if cfg.custom_css_files:
        config.custom_css_files.extend(cfg.custom_css_files)
    if cfg.module_order:
        config.module_order = cfg.module_order
    if cfg.extra_fn_clean_exts:
        config.fn_clean_exts = list(cfg.extra_fn_clean_exts) + config.fn_clean_exts
    if cfg.extra_fn_clean_trim:
        config.fn_clean_trim = list(cfg.extra_fn_clean_trim) + config.fn_clean_trim
    if cfg.preserve_module_raw_data is not None:
        config.preserve_module_raw_data = cfg.preserve_module_raw_data
    if cfg.data_dump_file_write_raw is not None:
        config.data_dump_file_write_raw = cfg.data_dump_file_write_raw
    if cfg.ai_summary is not None:
        config.ai_summary = cfg.ai_summary
    if cfg.ai_summary_full is not None:
        config.ai_summary = cfg.ai_summary_full
        config.ai_summary_full = cfg.ai_summary_full
    if cfg.ai_provider is not None:
        config.ai_provider = cfg.ai_provider
    if cfg.ai_model is not None:
        config.ai_model = cfg.ai_model
    if cfg.ai_custom_endpoint is not None:
        config.ai_custom_endpoint = cfg.ai_custom_endpoint
    if cfg.ai_custom_context_window is not None:
        config.ai_custom_context_window = cfg.ai_custom_context_window
    if cfg.ai_prompt_short is not None:
        config.ai_prompt_short = cfg.ai_prompt_short
    if cfg.ai_prompt_full is not None:
        config.ai_prompt_full = cfg.ai_prompt_full
    if cfg.no_ai is not None:
        config.no_ai = cfg.no_ai

    if config.development and "png" not in config.export_plot_formats:
        config.export_plot_formats.append("png")

    # Clean up analysis_dir if a string (interactive environment only)
    if analysis_dir:
        config.analysis_dir = list(analysis_dir)
    if cfg.file_list is not None:
        if len(config.analysis_dir) > 1:
            raise RunError("If --file-list is given, analysis_dir should have only one plain text file.")
        config.file_list = cfg.file_list

    if len(cfg.ignore) > 0:
        logger.debug(f"Ignoring files, directories and paths that match: {', '.join(cfg.ignore)}")
        config.fn_ignore_files.extend(cfg.ignore)
        config.fn_ignore_dirs.extend(cfg.ignore)
        config.fn_ignore_paths.extend(cfg.ignore)
    if len(cfg.ignore_samples) > 0:
        logger.debug(f"Ignoring sample names that match: {', '.join(cfg.ignore_samples)}")
        config.sample_names_ignore.extend(cfg.ignore_samples)
    if len(cfg.only_samples) > 0:
        logger.debug(f"Only including sample names that match: {', '.join(cfg.only_samples)}")
        config.sample_names_only_include.extend(cfg.only_samples)
    # Prep module configs
    report.top_modules = [m if isinstance(m, dict) else {m: {}} for m in config.top_modules]
    report.module_order = [m if isinstance(m, dict) else {m: {}} for m in config.module_order]

    if cfg.unknown_options:
        config.kwargs = cfg.unknown_options  # plug in command line options

    plugin_hooks.mqc_trigger("config_loaded")
    plugin_hooks.mqc_trigger("execution_start")
