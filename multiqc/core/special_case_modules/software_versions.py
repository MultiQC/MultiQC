"""Super Special-Case MultiQC module to produce report section on software versions"""

import logging
from textwrap import dedent
from typing import Dict, List

from multiqc import config
from multiqc import report
from multiqc.base_module import BaseMultiqcModule
from multiqc.types import Anchor
from multiqc.utils.material_icons import get_material_icon

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Software Versions",
            anchor=Anchor("multiqc_software_versions"),
            info="lists versions of software tools extracted from file contents.",
        )

        self.report_software_versions()
        self.write_software_versions_data_file()

    def report_software_versions(self):
        """Create section listing software versions."""
        content = self._make_versions_html(report.software_versions)
        self.add_section(name=None, content=content)

    @staticmethod
    def _make_versions_html(versions: Dict[str, Dict[str, List[str]]]) -> str:
        """Generate a tabular HTML output of all versions."""
        table_id = report.save_htmlid("mqc_versions_table")

        # Check if the Group column is identical to Software column
        groups_rows = []
        software_rows = []
        group_versions: Dict[str, List[str]]
        for group, group_versions in sorted(versions.items()):
            for tool, _ in sorted(group_versions.items()):
                groups_rows.append(group)
                software_rows.append(tool)
        ignore_groups = groups_rows == software_rows

        # Based on: https://github.com/nf-core/rnaseq/blob/3bec2331cac2b5ff88a1dc71a21fab6529b57a0f/modules/nf-core/custom/dumpsoftwareversions/templates/dumpsoftwareversions.py#L12
        header_rows = ["<th>Software</th>", "<th>Version</th>"]
        if not ignore_groups:
            header_rows.insert(0, f"<th>{config.versions_table_group_header}</th>")
        html = [
            dedent(
                f"""\
                <button type="button" class="mqc_table_copy_btn btn btn-outline-secondary btn-sm" data-clipboard-target="table#{table_id}">
                    {get_material_icon("mdi:content-copy", 16)} Copy table
                </button>
                <table class="table table-striped w-100 mqc_versions_table" id="{table_id}">
                    <thead>
                        <tr>{"".join(header_rows)}</tr>
                    </thead>
                """
            )
        ]
        for group, group_versions in sorted(versions.items()):
            html.append("<tbody>")
            tool_versions: List[str]
            for i, (tool, tool_versions) in enumerate(sorted(group_versions.items())):
                rows = [
                    f"<td>{tool}</td>",
                    f"<td><samp>{', '.join(list(map(str, tool_versions)))}</samp></td>",
                ]
                if not ignore_groups:
                    rows.insert(0, f"<td>{group if (i == 0) else ''}</td>")
                html.append(f"<tr>{''.join(rows)}</tr>")
            html.append("</tbody>")
        html.append("</table>")
        return "\n".join(html)

    @staticmethod
    def write_software_versions_data_file():
        """
        Write software versions to a file for downstream use
        """
        # Get rid of the default dicts and Version objects
        clean_software_versions: Dict[str, Dict[str, List[str]]] = {
            group: {software: list(map(str, svs)) for software, svs in versions.items()}
            for group, versions in report.software_versions.items()
        }

        # TSV only allows 2 levels of nesting.
        if config.data_format == "tsv":
            flat_software_versions: Dict[str, Dict[str, str]] = {
                group: {software: ", ".join(svs) for software, svs in versions.items()}
                for group, versions in clean_software_versions.items()
            }
            report.write_data_file(flat_software_versions, "multiqc_software_versions")

        else:
            report.write_data_file(clean_software_versions, "multiqc_software_versions")
