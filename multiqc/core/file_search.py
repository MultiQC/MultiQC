import logging
import os
from typing import Dict, List, Tuple

from multiqc.core.exceptions import RunError
from multiqc.utils import config, report

logger = logging.getLogger(__name__)


def file_search() -> Tuple[List, List]:
    """
    Search log files and set up the list of modules to run.
    """
    mod_keys = [list(m.keys())[0] for m in config.module_order]

    # Get the list of modules we want to run, in the order that we want them
    run_modules: List[Dict[str, Dict]] = [
        m for m in config.top_modules if list(m.keys())[0] in config.avail_modules.keys()
    ]
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
            raise RunError("No available modules to run!")
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
        raise RunError("No analysis modules specified!")
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
    if not required_logs_found(run_module_names):
        raise RunError()

    return run_modules, run_module_names


def required_logs_found(modules_with_logs):
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
