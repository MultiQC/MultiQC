import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.modules.ngsbits.mappingqc import parse_reports as mappingqc_parse_reports
from multiqc.modules.ngsbits.readqc import parse_reports as readqc_parse_reports


log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The ngs-bits module parses XML output generated for several tools in the ngs-bits collection:
    * [ReadQC](https://github.com/imgag/ngs-bits/blob/master/doc/tools/ReadQC.md) for statistics on FASTQ files,
    * [MappingQC](https://github.com/imgag/ngs-bits/blob/master/doc/tools/MappingQC.md) for statistics on BAM files
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="ngs-bits",
            anchor="ngsbits",
            href="https://github.com/imgag/ngs-bits",
            info="Calculating statistics from FASTQ, BAM, and VCF",
            doi="10.1093/bioinformatics/btx032",
        )

        # Call submodule functions
        n = dict()
        n["mappingqc"] = mappingqc_parse_reports(self)
        if n["mappingqc"] > 0:
            log.info(f"Found {n['mappingqc']} MappingQC reports")
        n["readqc"] = readqc_parse_reports(self)
        if n["readqc"] > 0:
            log.info(f"Found {n['readqc']} ReadQC reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
