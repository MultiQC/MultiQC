import logging
from typing import Dict, Set

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from . import family_sizes, grouping, consensus_families, yields, umis, consensus_stats, dedup, commands, review

_TOOL_MODULES = (
    family_sizes,
    grouping,
    consensus_families,
    yields,
    umis,
    consensus_stats,
    dedup,
    commands,
    review,
)

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    [fgumi](https://github.com/fulcrumgenomics/fgumi) is a fast Rust toolkit for working with Unique
    Molecular Identifiers (UMIs): extraction, grouping, consensus calling, correction and deduplication.
    """

    def __init__(self):
        super().__init__(
            name="fgumi",
            anchor="fgumi",
            href="https://github.com/fulcrumgenomics/fgumi",
            info="Fast Rust toolkit for UMI extraction, grouping, consensus calling, correction and deduplication.",
            # No DOI to cite // doi=
            license="MIT License",
            license_url="https://github.com/fulcrumgenomics/fgumi/blob/main/LICENSE",
        )

        self.samples_parsed_by_tool: Dict[str, Set[str]] = {}
        for mod in _TOOL_MODULES:
            tool_name = mod.__name__.rsplit(".", 1)[-1]
            samples = mod.parse_reports(self)
            self.samples_parsed_by_tool[tool_name] = samples
            if samples:
                log.info(f"Found {len(samples)} fgumi {tool_name} reports")

        if all(len(samples) == 0 for samples in self.samples_parsed_by_tool.values()):
            raise ModuleNoSamplesFound
