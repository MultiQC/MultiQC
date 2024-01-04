""" MultiQC module to parse output from vcftools """


import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .relatedness2 import Relatedness2Mixin
from .tstv_by_count import TsTvByCountMixin
from .tstv_by_qual import TsTvByQualMixin
from .tstv_summary import TsTvSummaryMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, Relatedness2Mixin, TsTvByCountMixin, TsTvByQualMixin, TsTvSummaryMixin):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="VCFTools",
            anchor="vcftools",
            href="https://vcftools.github.io",
            info="is a program for working with and reporting on VCF files.",
            doi="10.1093/bioinformatics/btr330",
        )

        n = dict()
        n["relatedness2"] = self.parse_relatedness2()

        n["tstv_by_count"] = self.parse_tstv_by_count()
        if n["tstv_by_count"] > 0:
            log.info(f"Found {n['tstv_by_count']} TsTv.count reports")

        n["tstv_by_qual"] = self.parse_tstv_by_qual()
        if n["tstv_by_qual"] > 0:
            log.info(f"Found {n['tstv_by_qual']} TsTv.qual reports")

        n["tstv_summary"] = self.parse_tstv_summary()
        if n["tstv_summary"] > 0:
            log.info(f"Found {n['tstv_summary']} TsTv.summary reports")

        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
