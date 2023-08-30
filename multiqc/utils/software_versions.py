#!/usr/bin/env python

""" Super Special-Case MultiQC module to produce report section on software versions """


import logging
import os

import yaml
import packaging.version

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import report

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Software Versions",
            anchor="multiqc_software_versions",
            info="lists versions of software tools extracted from file contents.",
        )

        self.report_software_versions()

    def report_software_versions(self):
        """Create section listing software versions."""
        content = "<dl class=dl-horizontal>\n"
        for tool_name in sorted(report.software_versions, key=lambda x: x.lower()):
            versions = [str(version) for version in report.software_versions[tool_name]]
            versions_string = "</code>, <code>".join(versions)
            content += f'  <dt style="text-align:left;">{tool_name}</dt><dd><code>{versions_string}</code></dd>\n'
        content += "</dl>\n"

        self.add_section(name=None, content=content)


def load_versions_from_config(config):
    """Try to load software versions from config"""

    def check_versions(versions):
        for version in versions:
            version, is_compliant = parse_version(str(version))
            if not is_compliant:
                log.debug(
                    f"'{software}' version '{version}' list in config.software_versions does not conform to PEP 440 format"
                )
            yield version

    software_versions = getattr(config, "software_versions", dict())

    file_name = config.version_fn_name
    if not os.path.isfile(file_name):
        file_name = file_name.replace(".yaml", ".yml")
    if os.path.isfile(file_name):
        with open(file_name) as f:
            try:
                log.debug("Reading software versions settings from: {}".format(file_name))
                software_versions.update(yaml.safe_load(f))
            except yaml.scanner.ScannerError as e:
                log.error("Error parsing versions YAML: {}".format(e))

    # Parse the versions in YAML and make sure that each software maps to a list of version strings
    for software in list(software_versions):
        versions = software_versions[software]
        if type(versions) != list:
            versions = [versions]

        # Check if versions are PEP 440 compliant and remove duplicates
        versions = list(set(check_versions(versions)))
        versions.sort(reverse=True)

        software_versions[software] = versions

    return software_versions


def find_matching_module(software_name: str, modules):
    """
    Find the module by name, ignoring case.
    """
    d = {m.name.lower(): m for m in modules}
    return d.get(software_name.lower())


def parse_version(version: str):
    """
    Check if version string is PEP 440 compliant to enable version normalization and proper ordering.
    Returns tuple with version and a boolean indicating if version is PEP 440 compliant.
    # - https://peps.python.org/pep-0440/
    """
    try:
        version = packaging.version.parse(version)
    except packaging.version.InvalidVersion:
        return version, False
    return version, True
