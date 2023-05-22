#!/usr/bin/env python

""" Super Special-Case MultiQC module to produce report section on software versions """


import logging

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
