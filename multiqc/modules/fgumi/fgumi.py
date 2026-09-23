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

    The module reads every metrics file fgumi writes. Files are recognized by their header row, so outputs
    with user-chosen names (for example `--stats` and `--metrics`) are found wherever they are:

    - `group`: family-size histogram, grouping metrics, position-group sizes
    - `simplex-metrics` / `duplex-metrics` (and `--metrics` on `simplex`, `duplex`, `codec`): family sizes,
      duplex AB x BA family sizes, yield by read pairs, UMI counts
    - `simplex` / `duplex` / `codec --stats`: consensus calling statistics
    - `correct`, `dedup` (metrics, family-size histogram, duplication ladder), `clip`, `filter`,
      `copy-umi`, `retag`, `downsample` histograms, and the `review` detail file

    The column contract for all of these is published in fgumi as `crates/fgumi-metrics/metric_columns.json`.

    #### Sample names

    Sample names come from the file name with fgumi's fixed suffixes (such as `.family_sizes.txt` or
    `.duplex_yield_metrics.txt`) removed, so all of one sample's `--metrics <prefix>` files share a name.

    #### fgbio GroupReadsByUmi histograms

    fgumi writes its family-size histogram with exactly the columns of fgbio GroupReadsByUmi's
    `--family-size-histogram`, so the two cannot be told apart. This module claims those files, so they
    appear under fgumi rather than fgbio.

    #### Large per-UMI files

    UMI counts, UMI correction and review files can exceed MultiQC's default file size limit, in which case
    they are skipped. Raise the limit to include them:

    ```yaml
    log_filesize_limit: 500000000
    ```
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
