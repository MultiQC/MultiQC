#!/usr/bin/env python

""" Super Special-Case MultiQC module to produce report section on software versions """


import logging
import os

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import report

import yaml

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Software Versions",
            anchor="multiqc_software_versions",
            info="Versions of software tools extracted from file contents.",
        )

        # Don't repeat the Custom Content name in the subtext
        if self.info or self.extra:
            self.intro = "<p>{}</p>{}".format(self.info, self.extra)

        self.report_software_versions()

    def report_software_versions(self):
        """Create section listing software versions."""
        content = "<dl class=dl-horizontal>\n"
        for tool_name in sorted(report.software_versions):
            versions = [str(version) for version in report.software_versions[tool_name]]
            versions_string = "</code>, <code>".join(versions)
            content += f"  <dt>{tool_name}</dt><dd><code>{versions_string}</code></dd>\n"
        content += "</dl>\n"

        self.add_section(name=None, content=content)


def load_versions_from_yaml(file_name):
    """Try to load software versions from any YAML file in current directory."""
    if not os.path.isfile(file_name):
        file_name = file_name.replace(".yaml", ".yml")

    if not os.path.isfile(file_name):
        return {}

    with open(file_name) as f:
        try:
            log.debug("Reading software versions settings from: {}".format(file_name))
            software_versions = yaml.safe_load(f)
        except yaml.scanner.ScannerError as e:
            log.error("Error parsing versions YAML: {}".format(e))
            return {}

    # Parse the versions in YAML and make sure that each software maps to a list of version strings
    for software in list(software_versions):
        versions = software_versions[software]
        if type(versions) != list:
            versions = [versions]
        versions = [str(version) for version in versions]
        software_versions[software] = versions

    return software_versions


def find_matching_module(software_name: str, modules):
    if software_name in modules:
        return software_name

    if software_name.capitalize() in modules:
        return software_name.capitalize()

    return None
