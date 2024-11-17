import glob
import logging
import os.path
from pathlib import Path
from typing import Dict, List, Tuple

from multiqc.core.exceptions import RunError, NoAnalysisFound
from multiqc import config, report

logger = logging.getLogger(__name__)


def file_search():
    """
    Search log files and set up the list of modules to run.
    """
    _make_analysis_file_list()

    mod_dicts_in_order, sp_keys = _module_list_to_search()
    report.search_files(sp_keys)

    return mod_dicts_in_order


def _make_analysis_file_list():
    """
    From config.file_list and config.analysis_dir, create a list of files
    to search in report.analysis_files
    """

    # Add files if --file-list option is given
    if config.file_list:
        file_list_path = Path(config.analysis_dir[0])
        with file_list_path.open() as in_handle:
            for line in in_handle:
                p = line.strip()
                if Path(p).exists():
                    report.analysis_files.append(str(Path(p).absolute()))
        if len(report.analysis_files) == 0:
            raise RunError(
                f"No files or directories were added from {file_list_path} using --file-list option."
                f"Please, check that {file_list_path} contains correct paths."
            )
    else:
        for path in config.analysis_dir:
            expanded_paths = list(glob.glob(str(path)))
            if not expanded_paths:
                logger.warning(f"No files found by path or pattern: '{path}'")
            for p in expanded_paths:  # Expand glob patterns
                report.analysis_files.append(str(p))

    if not report.analysis_files:
        raise NoAnalysisFound("No files found to analyse. Check that input files and directories exist")

    for p in report.analysis_files:
        logger.info(f"Search path: {os.path.abspath(p)}")


def include_or_exclude_modules(module_names: List[str]) -> List[str]:
    """
    Apply config.run_modules and config.exclude_modules filters
    """
    avail_modules = set(config.avail_modules.keys()) | {"software_versions"}
    if len(config.run_modules) > 0:
        unknown_modules = [m for m in config.run_modules if m not in avail_modules]
        if unknown_modules:
            logger.error(
                f"Module(s) in config.run_modules are unknown: {', '.join(unknown_modules)}. "
                f"Double check the case and dashes/underscores. Available modules IDs: {', '.join(config.avail_modules.keys())}"
            )
        if len(unknown_modules) == len(config.run_modules):
            raise RunError("No available modules to run!")
        module_names = [m for m in module_names if m in config.run_modules and m in avail_modules]
        if "custom_content" in config.run_modules:
            module_names.extend(config.custom_content_modules)
        logger.info(f"Only using modules: {', '.join(module_names)}")

    if len(config.exclude_modules) > 0:
        logger.info("Excluding modules '{}'".format("', '".join(config.exclude_modules)))
        if "general_stats" in config.exclude_modules:
            config.skip_generalstats = True
            config.exclude_modules = [x for x in config.exclude_modules if x != "general_stats"]
        module_names = [m for m in module_names if m not in config.exclude_modules]
    return module_names


def _module_list_to_search() -> Tuple[List[Dict[str, Dict]], List[str]]:
    """
    Get the list of modules we want to run, in the order that we want them.

    Return the list of search pattern keys, and the list of modules/entry points.
    """

    # Build initial list from report.module_order and config.top_modules
    mod_dicts_in_order: List[Dict[str, Dict]] = [
        m for m in report.top_modules if list(m.keys())[0] in config.avail_modules.keys()
    ]

    mod_order_keys = set(list(m.keys())[0] for m in report.module_order)
    mod_dicts_in_order.extend(
        [{m: {}} for m in config.avail_modules.keys() if m not in mod_order_keys and m not in mod_dicts_in_order]
    )
    mod_dicts_in_order.extend(
        [
            m
            for m in report.module_order
            if list(m.keys())[0] in config.avail_modules.keys()
            and list(m.keys())[0] not in [list(rm.keys())[0] for rm in mod_dicts_in_order]
        ]
    )

    # Always run software_versions module to collect version YAML files
    # Use config.skip_versions_section to exclude from report
    if "software_versions" not in mod_order_keys:
        mod_dicts_in_order.append({"software_versions": {}})

    mod_ids = include_or_exclude_modules([list(m.keys())[0] for m in mod_dicts_in_order])
    mod_dicts_in_order = [m for m in mod_dicts_in_order if list(m.keys())[0] in mod_ids]
    if len(mod_dicts_in_order) == 0:
        raise RunError("No analysis modules specified!")
    assert len(mod_dicts_in_order) == len(mod_ids)

    sp_keys = list(mod_ids)  # make a copy
    # Add custom content section names
    try:
        if "custom_content" in sp_keys:
            sp_keys = list(config.custom_data.keys()) + sp_keys
    except AttributeError:
        pass  # custom_data not in config

    # Always add software_versions
    if "software_versions" not in sp_keys:
        sp_keys.append("software_versions")

    logger.debug(f"Analysing modules: {', '.join(mod_ids)}")
    if sp_keys != mod_ids:
        logger.debug(f"Search keys: {', '.join(sorted(sp_keys))}")
    return mod_dicts_in_order, sp_keys
