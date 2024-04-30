import logging
import os
from typing import Dict, List

from multiqc.core.exceptions import RunError
from multiqc.utils import config, report

logger = logging.getLogger(__name__)


def _module_list_to_search() -> List[Dict[str, Dict]]:
    """
    Get the list of modules to search for files.
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
    if len(config.run_modules) > 0:
        unknown_modules = [m for m in config.run_modules if m not in config.avail_modules.keys()]
        if unknown_modules:
            logger.error(f"Module(s) in config.run_modules are unknown: {', '.join(unknown_modules)}")
        if len(unknown_modules) == len(config.run_modules):
            raise RunError("No available modules to run!")
        config.run_modules = [m for m in config.run_modules if m in config.avail_modules.keys()]
        run_modules = [m for m in run_modules if list(m.keys())[0] in config.run_modules]
        logger.info(f"Only using modules: {', '.join(config.run_modules)}")
    if len(config.exclude_modules) > 0:
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
        logger.info(f"Search path: {os.path.abspath(d)}")

    return run_modules


def file_search():
    """
    Search log files and set up the list of modules to run.
    """
    modules_to_search = _module_list_to_search()
    module_names = [list(m.keys())[0] for m in modules_to_search]

    report.search_files(module_names)

    return modules_to_search
