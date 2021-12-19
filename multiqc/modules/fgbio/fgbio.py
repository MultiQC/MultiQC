""" MultiQC module to parse output from fgbio """
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule

from .groupreadsbyumi import GroupReadsByUmiMixin
from . import ErrorRateByReadPosition


# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, GroupReadsByUmiMixin):
    """fgbio has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="fgbio",
            anchor="fgbio",
            target="fgbio",
            href="http://fulcrumgenomics.github.io/fgbio/",
            info=" is a command line toolkit for working with genomic and particularly next generation sequencing data..",
            # No publication / DOI // doi=
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()

        # GroupReadsByUmi
        n = dict()
        n["groupreadsbyumi"] = self.parse_groupreadsbyumi()
        if n["groupreadsbyumi"] > 0:
            log.info("Found {} groupreadsbyumi reports".format(n["groupreadsbyumi"]))

        # ErrorRateByReadPoosition
        n["errorratebyreadposition"] = ErrorRateByReadPosition.parse_reports(self)
        if n["errorratebyreadposition"] > 0:
            log.info("Found {} errorratebyreadposition reports".format(n["errorratebyreadposition"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
