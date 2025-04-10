"""Utility functions to handle software version reporting"""

import logging
import os
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import packaging.version
import yaml

from multiqc import config, report
from multiqc.types import ModuleId

# Initialise the logger
log = logging.getLogger(__name__)


def normalize_name(name: str):
    """Normalize software name to lowercase and remove spaces, underscores and dashes"""
    return name.lower().replace(" ", "").replace("-", "").replace("_", "")


def update_versions_from_config():
    """Update report with software versions from config if provided"""
    # Parse software version from config if provided
    versions_from_config: Dict[str, Dict[str, List[str]]] = load_versions_from_config()
    softwares: Dict[str, List[str]]
    for group, softwares in versions_from_config.items():
        # Try to find if the software is listed among the executed modules.
        # Unlisted software are still reported in the `Software Versions` section.
        module = find_matching_module(group, report.modules)

        # Map normalized module software names to the nicely formatted names.
        module_softwares = {}
        if module is not None:
            module_softwares = {normalize_name(m_software): m_software for m_software in module.versions.keys()}

            # Use the nicely formatted module name as group name
            group = module.name

        software_versions: List[str]
        for software_name, software_versions in softwares.items():
            # Update versions if the software is listed among the executed modules
            if module is not None and not config.disable_version_detection:
                # Update software name to match the module name format if found
                software_name = module_softwares.get(normalize_name(software_name), software_name)

                for version in software_versions:
                    module.add_software_version(str(version), software_name=software_name)

                # Get the updated software versions from the module
                software_versions = [version for _, version in module.versions[software_name]]

            # Add updated software versions to the report
            report.software_versions[group][software_name] = software_versions


def load_versions_from_config() -> Dict[str, Dict[str, List[str]]]:
    """Try to load software versions from config"""
    log.debug("Reading software versions from config.software_versions")
    versions_config: Dict[str, Dict[str, List[str]]]
    if not isinstance(config.software_versions, dict):
        log.error("Expected the `software_versions` config section to be a dictionary")
        versions_config = {}
    else:
        versions_config = validate_software_versions(config.software_versions)

    versions_from_files: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    for f in report.files.get(ModuleId("software_versions"), []):
        file_name = os.path.join(f["root"], f["fn"])
        with open(file_name) as fh:
            log.debug(f"Reading software versions settings from: {file_name}")
            try:
                versions_from_one_file = yaml.load(
                    fh,
                    # We need to be cautious when loading unquoted version strings from a YAML file.
                    # For instance, the version `1.10` will be parsed as a float by default, this converted
                    # into `1.1`. Passing yaml.BaseLoader explicitly makes YAML treat all scalar values
                    # as strings, so `1.10` will turn into a string `"1.10"` as we want.
                    # From https://pyyaml.org/wiki/PyYAMLDocumentation
                    #      BaseLoader(stream) does not resolve or support any tags
                    #      and constructs only basic Python objects: lists, dictionaries and Unicode strings.
                    Loader=yaml.BaseLoader,
                )
            except yaml.scanner.ScannerError as e:
                log.error(f"Error parsing versions YAML: {e}")

        if not isinstance(versions_from_one_file, dict):
            log.error(
                f"Expected the software versions file {file_name} to contain a dictionary structure, ignoring the file."
            )
            continue
        versions_from_one_file = validate_software_versions(versions_from_one_file)
        config.update_dict(versions_from_files, versions_from_one_file)

    # Aggregate versions listed in config and file
    config.update_dict(versions_config, versions_from_files)

    # Parse the aggregated versions
    for group in versions_config:
        softwares: Dict[str, List[str]] = versions_config[group]
        for tool, tool_versions in softwares.items():
            # Try and convert version to packaging.versions.Version object and remove duplicates
            ver_and_verstr = list(set((parse_version(version), version) for version in tool_versions))
            ver_and_verstr = sort_versions(ver_and_verstr)
            softwares[tool] = [v for _, v in ver_and_verstr]

    return versions_config


def validate_software_versions(versions_config: Dict[str, Any]) -> Dict[str, Dict[str, List[str]]]:
    """
    Validate software versions input from config file

    Two main formats are supported:

    (1) A dict of dicts, where the first level is the group name and the second level is the software name.
        The software name maps to a list of versions or a single version string.

        group_name:
            software_name:
                - version1
                - version2
        other_group:
            other_software: version3

    (2) A dict of lists, where the key is the software name and the value is a list of versions or a single version string.
        In this case the group name is set to the software name.

        software_name:
            - version1
            - version2
        other_software: version3

    It is also possible to mix the two formats, but this is not recommended.

    Returns a dict of dicts of list in the format (1).
    """

    def _filter_list(lst: List[str]) -> List[str]:
        """Remove all non-string version tags"""
        fixed_lst = []
        for item in lst:
            if not isinstance(item, str):
                log.error(
                    f"Version must be a string, got '{type(item).__name__}': '{item}'. Consider wrapping the value in quotes: '\"{item}\"'"
                )
            else:
                fixed_lst.append(item)
        return fixed_lst

    output: Dict[str, Dict[str, List]] = defaultdict(lambda: defaultdict(list))

    for level1_key, level1_values in versions_config.items():
        group = level1_key
        software = level1_key

        # Check if the input is in format (1)
        if isinstance(level1_values, dict):
            for level2_key, versions in level1_values.items():
                software = level2_key
                if not isinstance(versions, list):
                    versions = [versions]
                versions = _filter_list(versions)
                if versions:
                    output[group][software] = versions

        # Check if the input is in format (2)
        else:
            versions = level1_values
            if not isinstance(versions, list):
                versions = [versions]
            versions = _filter_list(versions)
            if versions:
                output[group][software] = versions

    return output


def sort_versions(
    ver_and_verstr: Sequence[Tuple[Optional[packaging.version.Version], str]],
) -> List[Tuple[Optional[packaging.version.Version], str]]:
    """
    Sort list of versions in descending order. Accepts list with both strings and packaging.version.Version
    objects.
    """
    version_parsed = [(vobj, vstr) for vobj, vstr in ver_and_verstr if vobj is not None]
    version_not_parsed = [(vobj, vstr) for vobj, vstr in ver_and_verstr if vobj is None]
    versions = sorted(version_parsed) + sorted(version_not_parsed)
    return list(versions)


def find_matching_module(software_name: str, modules):
    """
    Find the module by name
    """
    d = {normalize_name(m.name): m for m in modules}
    return d.get(normalize_name(software_name))


def parse_version(version: str) -> Optional[packaging.version.Version]:
    """
    Check if version string is PEP 440 compliant to enable version normalization and proper ordering.
    Returns tuple with version and a boolean indicating if version is PEP 440 compliant.
    # - https://peps.python.org/pep-0440/
    """
    try:
        return packaging.version.parse(version)
    except packaging.version.InvalidVersion:
        return None
