import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from . import alignment, basic, gcbias, hybcap, isize, wgs

TOOLS = [
    m.__name__.split(".")[-1]
    for m in (
        alignment,
        basic,
        gcbias,
        hybcap,
        isize,
        wgs,
    )
]

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    [Riker](https://github.com/fulcrumgenomics/riker) is a fast Rust toolkit for
    sequencing QC metrics that ports many of the most widely-used tools from Picard
    with cleaner output and better performance.

    Supported subtools:

    - `alignment` (equivalent to Picard's `CollectAlignmentSummaryMetrics`)
    - `basic` (equivalent to `CollectBaseDistributionByCycle`, `MeanQualityByCycle`,
      and `QualityScoreDistribution`)
    - `gcbias` (equivalent to `CollectGcBiasMetrics`)
    - `hybcap` (equivalent to `CollectHsMetrics`)
    - `isize` (equivalent to `CollectInsertSizeMetrics`)
    - `wgs` (equivalent to `CollectWgsMetrics`)

    Riker emits plain TSV files with `sample` as the first column and snake_case
    column names; this module parses those files directly. Per-target and per-base
    coverage outputs from `hybcap` (only emitted with `--per-target-coverage`)
    are not parsed in this version. The `error` subtool is also not yet supported.

    #### Coverage thresholds

    Riker emits cumulative coverage fractions for `wgs` (e.g. `frac_bases_at_30x`)
    and `hybcap` (e.g. `frac_target_bases_30x`). The threshold shown in the General
    Statistics table can be customised:

    ```yaml
    riker_config:
      general_stats_target_coverage:
        - 10
        - 30
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Riker",
            anchor="riker",
            href="https://github.com/fulcrumgenomics/riker",
            info="Fast Rust toolkit that ports key sequencing QC tools from Picard.",
            # No DOI to cite // doi=
        )

        self.samples_parsed_by_tool = dict()

        for tool_name in TOOLS:
            log.debug(f"Running riker tool {tool_name}")
            mod = globals()[tool_name]
            samples = mod.parse_reports(self)
            self.samples_parsed_by_tool[tool_name] = samples
            if len(samples) > 0:
                log.info(f"Found {len(samples)} riker {tool_name} reports")

        if all(len(v) == 0 for v in self.samples_parsed_by_tool.values()):
            raise ModuleNoSamplesFound
