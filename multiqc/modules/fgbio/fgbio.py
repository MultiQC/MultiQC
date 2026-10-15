import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.modules.fgbio.collect_duplex_seq_metrics import collect_duplex_seq_metrics
from multiqc.modules.fgbio.error_rate_by_read_position import error_rate_by_read_position
from multiqc.modules.fgbio.group_reads_by_umi import run_group_reads_by_umi

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module currently supports tool outputs from:

    - [GroupReadsByUmi](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html)

      fgumi writes this family-size histogram with identical columns, and the fgumi module reads it too. The fgbio
      module keeps reporting it unless its directory also holds files only fgumi writes, in which case it appears
      under fgumi instead. Set `fgumi_config: {family_sizes_module: fgbio}` to always keep it here (see the fgumi
      module docs).
    - [ErrorRateByReadPosition](http://fulcrumgenomics.github.io/fgbio/tools/latest/ErrorRateByReadPosition.html)
    - [CollectDuplexSeqMetrics](http://fulcrumgenomics.github.io/fgbio/tools/latest/CollectDuplexSeqMetrics.html)
    """

    def __init__(self):
        super().__init__(
            name="fgbio",
            anchor="fgbio",
            target="fgbio",
            href="http://fulcrumgenomics.github.io/fgbio/",
            info="Processing and evaluating data containing UMIs",
            # No publication / DOI // doi=
            license="MIT License",
            license_url="https://github.com/fulcrumgenomics/fgbio/blob/master/LICENSE",
        )

        n = {}

        # GroupReadsByUmi
        n["groupreadsbyumi"] = run_group_reads_by_umi(self)
        if n["groupreadsbyumi"] > 0:
            log.info(f"Found {n['groupreadsbyumi']} groupreadsbyumi reports")

        # ErrorRateByReadPosition
        n["errorratebyreadposition"] = error_rate_by_read_position(self)
        if n["errorratebyreadposition"] > 0:
            log.info(f"Found {n['errorratebyreadposition']} errorratebyreadposition reports")

        # CollectDuplexSeqMetrics
        n["collectduplexseqmetrics"] = collect_duplex_seq_metrics(self)
        if n["collectduplexseqmetrics"] > 0:
            log.info(f"Found {n['collectduplexseqmetrics']} collectduplexseqmetrics reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
