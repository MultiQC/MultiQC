""" Utility functions to handle software version reporting """

import logging
import os
from collections import defaultdict
from typing import List, Dict

import packaging.version
import yaml

from multiqc.utils import report as mqc_report

# Initialise the logger
log = logging.getLogger(__name__)


def normalize_name(name: str):
    """Normalize software name to lowercase and remove spaces, underscores and dashes"""
    return name.lower().replace(" ", "").replace("-", "").replace("_", "")


def update_versions_from_config(config, report):
    """Update report with software versions from config if provided"""
    # Parse software version from config if provided
    versions_from_config = load_versions_from_config(config)
    for group, softwares in versions_from_config.items():
        # Try to find if the software is listed among the executed modules.
        # Unlisted software are still reported in the `Software Versions` section.
        module = find_matching_module(group, report.modules_output)

        # Map normalized module software names to the nicely formatted names.
        module_softwares = {}
        if module is not None:
            module_softwares = {normalize_name(m_software): m_software for m_software in module.versions}

            # Use the nicely formatted module name as group name
            group = module.name

        for software, versions in softwares.items():
            # Update versions if the software is listed among the executed modules
            if module is not None and not config.disable_version_detection:
                # Update software name to match the module name format if found
                software = module_softwares.get(normalize_name(software), software)

                for version in versions:
                    module.add_software_version(str(version), software_name=software)

                # Get the updated software versions from the module
                versions = module.versions[software]

            # Add updated software versions to the report
            report.software_versions[group][software] = versions


def load_versions_from_config(config):
    """Try to load software versions from config"""
    log.debug("Reading software versions from config.software_versions")
    versions_config = getattr(config, "software_versions", defaultdict(lambda: defaultdict(list)))
    if not isinstance(versions_config, dict):
        log.error("Expected the `software_versions` config section to be a dictionary")
        versions_config = {}
    else:
        versions_config = validate_software_versions(versions_config)

    versions_from_files = defaultdict(lambda: defaultdict(list))
    for f in mqc_report.files.get("software_versions", []):
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
                f"Expected the software versions file {file_name} to contain a dictionary structure, "
                f"ignoring the file."
            )
            continue
        versions_from_one_file = validate_software_versions(versions_from_one_file)
        config.update_dict(versions_from_files, versions_from_one_file)

    # Aggregate versions listed in config and file
    config.update_dict(versions_config, versions_from_files)

    # Parse the aggregated versions
    for group in versions_config:
        softwares = versions_config[group]
        for tool in softwares:
            # Try and convert version to packaging.versions.Version object and remove duplicates
            versions = list(set([parse_version(version) for version in softwares[tool]]))
            versions = sort_versions(versions)
            softwares[tool] = versions

    return versions_config


def validate_software_versions(versions_config: Dict) -> Dict[str, Dict]:
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

    output = defaultdict(lambda: defaultdict(list))

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


def sort_versions(versions):
    """
    Sort list of versions in descending order. Accepts list with both strings and packaging.version.Version
    objects.
    """
    version_objs = [v for v in versions if isinstance(v, packaging.version.Version)]
    version_strs = [v for v in versions if not isinstance(v, packaging.version.Version)]
    versions = sorted(version_objs) + sorted(version_strs)
    return versions


def find_matching_module(software_name: str, modules):
    """
    Find the module by name
    """
    d = {normalize_name(m.name): m for m in modules}
    return d.get(normalize_name(software_name))


def parse_version(version: str):
    """
    Check if version string is PEP 440 compliant to enable version normalization and proper ordering.
    Returns tuple with version and a boolean indicating if version is PEP 440 compliant.
    # - https://peps.python.org/pep-0440/
    """
    try:
        return packaging.version.parse(version)
    except packaging.version.InvalidVersion:
        return str(version)
