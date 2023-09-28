#!/usr/bin/env python

""" Super Special-Case MultiQC module to produce report section on software versions """


import logging
import os
from collections import defaultdict

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
    software_versions_config = getattr(config, "software_versions", defaultdict(lambda: defaultdict(list)))
    software_versions_config, is_valid = validate_software_versions(software_versions_config)

    if not is_valid:
        log.error("Software versions loaded config.software_versions is not in a valid format")

    software_versions_file = defaultdict(lambda: defaultdict(list))
    for f in mqc_report.files.get("software_versions", []):
        file_name = os.path.join(f["root"], f["fn"])
        with open(file_name) as f:
            try:
                log.debug("Reading software versions settings from: {}".format(file_name))
                software_versions_file_tmp = yaml.safe_load(f)
            except yaml.scanner.ScannerError as e:
                log.error("Error parsing versions YAML: {}".format(e))

        software_versions_file_tmp, is_valid = validate_software_versions(software_versions_file_tmp)

        if not is_valid:
            log.error("Software versions loaded from {} is not in a valid format".format(file_name))
            continue

        software_versions_file = merge(software_versions_file, software_versions_file_tmp)

    # Aggregate versions listed in config and file
    software_versions = merge(software_versions_config, software_versions_file)

    # Parse the aggregated versions
    for group in list(software_versions):
        softwares = software_versions[group]
        for tool in list(softwares):
            versions = softwares[tool]

            # Try and convert version to packaging.versions.Version object and remove duplicates
            versions = list(set([parse_version(version) for version in versions]))
            versions = sort_versions(versions)

            softwares[tool] = versions

    return software_versions


def validate_software_versions(input):
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

    Returns a dict of dicts of list in the format (1) and a boolean indicating if the input is valid.
    """

    def all_strings(lst):
        return all(isinstance(item, str) for item in lst)

    output = defaultdict(lambda: defaultdict(list))
    if not isinstance(input, dict):
        return output, False

    for level1_key, level1_values in input.items():
        group = level1_key
        software = level1_key

        # Check if the input is in format (1)
        if isinstance(level1_values, dict):
            for level2_key, versions in level1_values.items():
                software = level2_key
                if isinstance(versions, str):
                    versions = [versions]

                if not isinstance(versions, list):
                    return output, False

                if not all_strings(versions):
                    return output, False

                output[group][software] = versions

        # Check if the input is in format (2)
        elif isinstance(level1_values, (list, str)):
            versions = level1_values
            if isinstance(versions, str):
                versions = [versions]

            if not all_strings(versions):
                return output, False

            output[group][software] = versions
        else:
            return output, False

    return output, True


def merge(a: dict, b: dict, path=None):
    """Merge two dict of dicts recursively"""
    # source: https://stackoverflow.com/a/7205107
    if path is None:
        path = []

    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif isinstance(a[key], list) and isinstance(b[key], list):
                a[key].extend(b[key])
            else:
                raise Exception("Conflict at " + ".".join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a


def sort_versions(versions):
    """Sort list of versions in descending order. Accepts list with both strings and packaging.version.Version objects."""
    try:
        versions.sort()
    except TypeError:
        # If there is a mix, sort all as strings
        versions.sort(key=str)
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
