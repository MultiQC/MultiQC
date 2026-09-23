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

    The column contract for all of these is fgumi's `crates/fgumi-metrics/metric_columns.json`. The `copy-umi`
    and `retag` metrics and the headered `filter --stats` layout need an fgumi release after 0.7.0; the headerless
    `filter --stats` layout written by fgumi 0.7.0 and earlier is also read.

    #### Sample names

    Sample names come from the file name with fgumi's fixed suffixes (such as `.family_sizes.txt` or
    `.duplex_yield_metrics.txt`) removed, so all of one sample's `--metrics <prefix>` files share a name.

    #### fgbio outputs

    Several fgumi metrics files have exactly the columns of their fgbio equivalents, so the two cannot be
    told apart, and this module also reads (and reports under fgumi) the fgbio versions of:

    - GroupReadsByUmi: `--family-size-histogram` and `--grouping-metrics` (the family-size histogram was
      previously shown by the fgbio module; it now appears under fgumi instead)
    - CollectDuplexSeqMetrics: family sizes, duplex family sizes, duplex yield, UMI counts, duplex UMI counts
    - CorrectUmis: `--metrics`
    - ClipBam: `--metrics`

    #### Duplex heatmaps

    Each sample gets AB x BA family-size heatmaps unless there are more than 10 duplex samples, in which case
    only the cross-sample plot is drawn. Change the limit with:

    ```yaml
    fgumi_config:
      max_duplex_heatmap_samples: 10
    ```

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
