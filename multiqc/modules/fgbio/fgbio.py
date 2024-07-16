"""MultiQC module to parse output from fgbio"""

import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from . import ErrorRateByReadPosition
from .groupreadsbyumi import GroupReadsByUmiMixin

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, GroupReadsByUmiMixin):
    """
    The fgbio MultiQC module currently supports tool the following outputs:

    - [GroupReadsByUmi](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html)
    - [ErrorRateByReadPosition](http://fulcrumgenomics.github.io/fgbio/tools/latest/ErrorRateByReadPosition.html)
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="fgbio",
            anchor="fgbio",
            target="fgbio",
            href="http://fulcrumgenomics.github.io/fgbio/",
            info="Processing and evaluating data containing UMIs",
            # No publication / DOI // doi=
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()

        # GroupReadsByUmi
        n = dict()
        n["groupreadsbyumi"] = self.parse_groupreadsbyumi()
        if n["groupreadsbyumi"] > 0:
            log.info(f"Found {n['groupreadsbyumi']} groupreadsbyumi reports")

        # ErrorRateByReadPoosition
        n["errorratebyreadposition"] = ErrorRateByReadPosition.parse_reports(self)
        if n["errorratebyreadposition"] > 0:
            log.info(f"Found {n['errorratebyreadposition']} errorratebyreadposition reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
