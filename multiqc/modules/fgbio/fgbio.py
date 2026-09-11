import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from multiqc.modules.fgbio.clip_bam import run_clip_bam
from multiqc.modules.fgbio.error_rate_by_read_position import error_rate_by_read_position
from multiqc.modules.fgbio.group_reads_by_umi import run_group_reads_by_umi

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module currently supports tool the following outputs:

    - [GroupReadsByUmi](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html)
    - [ErrorRateByReadPosition](http://fulcrumgenomics.github.io/fgbio/tools/latest/ErrorRateByReadPosition.html)
    - [ClipBam](http://fulcrumgenomics.github.io/fgbio/tools/latest/ClipBam.html)

    For `ClipBam`, the module reads the metrics file written with `--metrics` and reports, per
    read type, how many reads and bases were clipped and for which reason. Note that `bases` in
    that file counts the aligned bases left after clipping, so the percentages of bases clipped
    are computed against `bases + bases_clipped_post`, the bases present before ClipBam ran.
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

        # GroupReadsByUmi
        n = dict()
        n["groupreadsbyumi"] = run_group_reads_by_umi(self)
        if n["groupreadsbyumi"] > 0:
            log.info(f"Found {n['groupreadsbyumi']} groupreadsbyumi reports")

        # ErrorRateByReadPoosition
        n["errorratebyreadposition"] = error_rate_by_read_position(self)
        if n["errorratebyreadposition"] > 0:
            log.info(f"Found {n['errorratebyreadposition']} errorratebyreadposition reports")

        # ClipBam
        n["clipbam"] = run_clip_bam(self)
        if n["clipbam"] > 0:
            log.info(f"Found {n['clipbam']} clipbam reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
