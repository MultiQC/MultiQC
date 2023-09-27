#!/usr/bin/env python

""" Super Special-Case MultiQC module to produce report section on software versions """


import logging
from textwrap import dedent

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import config as mqc_config
from multiqc.utils import report as mqc_report
from multiqc.utils import util_functions

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
        self.write_software_versions_data_file()

    def report_software_versions(self):
        """Create section listing software versions."""
        content = self._make_versions_html(mqc_report.software_versions)
        self.add_section(name=None, content=content)

    @staticmethod
    def _make_versions_html(versions):
        """Generate a tabular HTML output of all versions."""
        table_id = mqc_report.save_htmlid("mqc_versions_table")

        # Check if the Group column is identical to Software column
        groups_rows = []
        software_rows = []
        for group, tmp_versions in sorted(versions.items()):
            for tool, _ in sorted(tmp_versions.items()):
                groups_rows.append(group)
                software_rows.append(tool)
        ignore_groups = groups_rows == software_rows

        # Based on: https://github.com/nf-core/rnaseq/blob/3bec2331cac2b5ff88a1dc71a21fab6529b57a0f/modules/nf-core/custom/dumpsoftwareversions/templates/dumpsoftwareversions.py#L12
        header_rows = ["<th>Software</th>", "<th>Version</th>"]
        if not ignore_groups:
            header_rows.insert(0, f"<th>{mqc_config.versions_table_group_header}</th>")
        html = [
            dedent(
                f"""\
                <button type="button" class="mqc_table_copy_btn btn btn-default btn-sm" data-clipboard-target="#{table_id}">
                    <span class="glyphicon glyphicon-copy"></span> Copy table
                </button>
                <table class="table mqc_versions_table" id="{table_id}">
                    <thead>
                        <tr>{''.join(header_rows)}</tr>
                    </thead>
                """
            )
        ]
        for group, tmp_versions in sorted(versions.items()):
            html.append("<tbody>")
            for i, (tool, versions) in enumerate(sorted(tmp_versions.items())):
                rows = [
                    f"<td>{tool}</td>",
                    f"<td><samp>{', '.join(list(map(str, versions)))}</samp></td>",
                ]
                if not ignore_groups:
                    rows.insert(0, f"<td>{group if (i == 0) else ''}</td>")
                html.append(f"<tr>{''.join(rows)}</tr>")
            html.append("</tbody>")
        html.append("</table>")
        return "\n".join(html)

    @staticmethod
    def write_software_versions_data_file():
        """Write software versions to a file for downstream use."""
        # Get rid of the default dicts and Version objects
        flat_software_versions = {
            group: {software: list(map(str, software_versions)) for software, software_versions in versions.items()}
            for group, versions in mqc_report.software_versions.items()
        }
        # TSV only allows 2 levels of nesting.
        if mqc_config.data_format == "tsv":
            flat_software_versions = {
                group: {software: ", ".join(software_versions) for software, software_versions in versions.items()}
                for group, versions in flat_software_versions.items()
            }
        # Write to a file for downstream use
        util_functions.write_data_file(flat_software_versions, "multiqc_software_versions")
