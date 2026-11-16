"""MultiQC module to render structured tool citations."""

import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .citation import Citation, render_bibliography, render_inline
from .csl import parse_csl

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Renders structured software citations for the tools used in an analysis.

    This module does not parse a tool's output. Instead it consumes a CSL-JSON
    citations file (`*.csl.json`) written by a workflow or pipeline (for example
    by nf-core pipelines via nf-core-utils), and renders two sections: a methods
    sentence listing each tool with the version that ran and a short citation,
    and a full bibliography.

    The file is a JSON array of
    [CSL-JSON](https://citeproc-js.readthedocs.io/en/latest/csl-json/) items. A
    `custom` object on each item carries the tool name and the version that
    actually ran (`{"custom": {"tool": "fastqc", "version": "0.12.1"}}`).

    Short citations use the Harvard short form ("Surname (year)", or
    "Surname et al. (year)" for several authors) and link to the DOI when
    present, otherwise the tool homepage. Tools with only a DOI and no
    publication are identified by tool name rather than rendered anonymously.
    """

    def __init__(self):
        # No doi= here: this module renders citations for many tools, it is not
        # itself a wrapper around one publication.
        super().__init__(
            name="Citations",
            anchor="citations",
            info="Software citations for the tools used in this analysis.",
        )

        citations = self._parse_citations()
        if not citations:
            raise ModuleNoSamplesFound

        # Required by linting even though this module has no own version.
        self.add_software_version(None)

        self.add_section(
            name="Methods",
            anchor="citations-methods",
            description="Tools run in this workflow, with the version used and a short citation.",
            content=render_inline(citations),
        )
        self.add_section(
            name="Bibliography",
            anchor="citations-bibliography",
            description="Full references for the tools cited above.",
            content=render_bibliography(citations),
        )

        self.write_data_file({c.tool: c.to_dict() for c in citations}, "multiqc_citations")

    def _parse_citations(self):
        """Collect citations across all matched files, deduped by tool.

        Tool names act as sample names, so `--ignore-samples` patterns filter
        tools out of the citation list.
        """
        by_tool: Dict[str, Citation] = {}
        for f in self.find_log_files("citations"):
            for citation in parse_csl(f["f"], f["fn"]):
                by_tool.setdefault(citation.tool, citation)
            self.add_data_source(f)

        by_tool = self.ignore_samples(by_tool)

        # Deterministic order for stable reports.
        return [by_tool[tool] for tool in sorted(by_tool, key=str.lower)]
