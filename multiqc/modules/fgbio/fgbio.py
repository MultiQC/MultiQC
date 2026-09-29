import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from multiqc.modules.fgbio.clip_bam import run_clip_bam
from multiqc.modules.fgbio.error_rate_by_read_position import error_rate_by_read_position
from multiqc.modules.fgbio.group_reads_by_umi import run_group_reads_by_umi
from multiqc.modules.fgumi.util import found_by, shared_file_evidence, shared_file_module

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

    fgumi writes the GroupReadsByUmi family-size histogram and the ClipBam metrics with identical
    columns, and the fgumi module reads them too. The fgbio module keeps reporting them unless their
    directory also holds files only fgumi writes, in which case they appear under fgumi instead. Set
    `fgumi_config: {shared_files_module: fgbio}` to always keep them here (see the fgumi module docs).
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
        evidence = shared_file_evidence()
        n["clipbam"] = len(
            run_clip_bam(
                self,
                "fgbio/clipbam",
                skip=lambda f: found_by(f, "fgumi/clip") and shared_file_module(f, evidence) == "fgumi",
            )
        )
        if n["clipbam"] > 0:
            log.info(f"Found {n['clipbam']} clipbam reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
