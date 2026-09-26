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

    fgumi publishes the columns of every one of these files in `crates/fgumi-metrics/metric_columns.json`, starting
    with the first release after 0.7.0. That release is also the first to write `copy-umi` and `retag` metrics and a
    headered `filter --stats` file. The headerless `filter --stats` file from fgumi 0.7.0 and earlier is read too.

    #### Sample names

    Sample names come from the file name with fgumi's fixed suffixes (such as `.family_sizes.txt` or
    `.duplex_yield_metrics.txt`) removed, so all of one sample's `--metrics <prefix>` files share a name. The stage
    that `fgumi runall --all-metrics <prefix>` adds before the suffix (as in `<prefix>.duplex.umi_counts.txt`) is
    removed too. With `--fullnames`, nothing is removed.

    #### fgbio outputs

    Several fgumi metrics files have exactly the columns of their fgbio equivalents. The two can't be told apart,
    so this module also reads the fgbio versions and reports them under fgumi:

    - GroupReadsByUmi: `--family-size-histogram` and `--grouping-metrics` (for the family-size histogram, see
      below)
    - CollectDuplexSeqMetrics: family sizes, duplex family sizes, duplex yield, UMI counts, duplex UMI counts
    - CorrectUmis: `--metrics`
    - ClipBam: `--metrics`
    - ReviewConsensusVariants: the `<output>.txt` detail file

    Of these, only the family-size histogram is also read by the fgbio module, and only one of the two modules
    reports each file. The histogram's own directory decides, or the nearest parent directory if its own has no
    evidence: fgumi when that directory holds files only fgumi writes (for example position group sizes, dedup,
    filter, copy-umi, retag or downsample metrics) and no files only fgbio writes (ErrorRateByReadPosition), and
    fgbio otherwise. With no evidence anywhere, fgbio reports it. Choose the module for every file with:

    ```yaml
    fgumi_config:
      family_sizes_module: fgumi # or fgbio
    ```

    #### Duplex heatmaps

    Each sample gets an AB x BA family-size heatmap of all families, plus one of only the families with reads on
    both strands when some families have no BA reads (otherwise the two would be identical). With more than 5
    duplex samples only the cross-sample plot is drawn. Change the limit with:

    ```yaml
    fgumi_config:
      max_duplex_heatmap_samples: 5
    ```

    #### Large per-UMI files

    UMI counts, UMI correction and review files can exceed MultiQC's default file size limit, in which case
    they are skipped. Raise the limit to include them (summarizing a large review file keeps every distinct
    consensus read name in memory, roughly 100 bytes per read):

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
