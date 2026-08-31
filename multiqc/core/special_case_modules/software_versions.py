"""Super Special-Case MultiQC module to produce report section on software versions"""

import logging
from textwrap import dedent
from typing import Dict, List, Optional

from multiqc import config
from multiqc import report
from multiqc.base_module import BaseMultiqcModule
from multiqc.types import Anchor, SoftwareVersionMetadata
from multiqc.utils.material_icons import get_material_icon

# Initialise the logger
log = logging.getLogger(__name__)


def _license_html(meta: Optional[SoftwareVersionMetadata]) -> str:
    """Render the License table cell, linking the license name to its URL if available."""
    if meta is None or (not meta.license and not meta.license_url):
        return ""
    label = meta.license or meta.license_url
    if meta.license_url:
        return f'<a href="{meta.license_url}" target="_blank">{label}</a>'
    return str(label)


def _doi_html(meta: Optional[SoftwareVersionMetadata]) -> str:
    """Render the DOI table cell as one or more links to doi.org."""
    if meta is None or not meta.doi:
        return ""
    links = [f'<a href="https://doi.org/{doi}" target="_blank">{doi}</a>' for doi in meta.doi]
    return "; ".join(links)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="Software Versions",
            anchor=Anchor("multiqc_software_versions"),
            info="lists versions of software tools extracted from file contents.",
            # Internal MultiQC section, not an external tool: no license to declare
            license=None,
            license_url=None,
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
        metadata = report.software_versions_metadata

        # Check if the Group column is identical to Software column
        groups_rows = []
        software_rows = []
        group_versions: Dict[str, List[str]]
        for group, group_versions in sorted(versions.items()):
            for tool, _ in sorted(group_versions.items()):
                groups_rows.append(group)
                software_rows.append(tool)
        ignore_groups = groups_rows == software_rows

        # Only show the optional License / DOI columns if any software provides them
        show_license = any(
            (meta.license or meta.license_url)
            for group, group_versions in versions.items()
            for tool in group_versions
            if (meta := metadata.get(group, {}).get(tool)) is not None
        )
        show_doi = any(
            meta.doi
            for group, group_versions in versions.items()
            for tool in group_versions
            if (meta := metadata.get(group, {}).get(tool)) is not None
        )

        # Based on: https://github.com/nf-core/rnaseq/blob/3bec2331cac2b5ff88a1dc71a21fab6529b57a0f/modules/nf-core/custom/dumpsoftwareversions/templates/dumpsoftwareversions.py#L12
        header_rows = ["<th>Software</th>", "<th>Version</th>"]
        if not ignore_groups:
            header_rows.insert(0, f"<th>{config.versions_table_group_header}</th>")
        if show_license:
            header_rows.append("<th>License</th>")
        if show_doi:
            header_rows.append("<th>DOI</th>")
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
                meta = metadata.get(group, {}).get(tool)
                rows = [
                    f"<td>{tool}</td>",
                    f"<td><samp>{', '.join(list(map(str, tool_versions)))}</samp></td>",
                ]
                if not ignore_groups:
                    rows.insert(0, f"<td>{group if (i == 0) else ''}</td>")
                if show_license:
                    rows.append(f"<td>{_license_html(meta)}</td>")
                if show_doi:
                    rows.append(f"<td>{_doi_html(meta)}</td>")
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
