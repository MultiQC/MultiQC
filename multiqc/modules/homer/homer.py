import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .findpeaks import FindPeaksReportMixin
from .tagdirectory import TagDirReportMixin

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, FindPeaksReportMixin, TagDirReportMixin):
    """
    HOMER contains many useful tools for analyzing ChIP-Seq, GRO-Seq, RNA-Seq, DNase-Seq, Hi-C and numerous
    other types of functional genomics sequencing data sets. The module currently only parses output from the
    `findPeaks` and `TagDirectory` tools. If you would like support to be added for other HOMER tools please
    open a [new issue](https://github.com/MultiQC/MultiQC/issues/new) on the MultiQC GitHub page.

    #### FindPeaks

    The HOMER findPeaks MultiQC module parses the summary statistics found at the top
    of HOMER peak files. Three key statistics are shown in the General Statistics table,
    all others are saved to `multiqc_data/multiqc_homer_findpeaks.txt`.

    #### TagDirectory

    The HOMER tag directory submodule parses output from files
    [tag directory](http://homer.ucsd.edu/homer/ngs/tagDir.html) output files, generating
    a number of diagnostic plots.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="HOMER",
            anchor="homer",
            href="http://homer.ucsd.edu/homer/",
            info="Motif discovery and next-gen sequencing analysis.",
            doi="10.1016/j.molcel.2010.05.004",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # Set up data structures
        self.tagdir_data = {
            "GCcontent": {},
            "restriction": {},
            "restriction_norm": {},
            "length": {},
            "taginfo_total": {},
            "taginfo_total_norm": {},
            "taginfo_uniq": {},
            "taginfo_uniq_norm": {},
            "FreqDistribution": {},
            "header": {},
            "interChr": {},
        }
        # Call submodule functions

        n["findpeaks"] = self.parse_homer_findpeaks()
        if n["findpeaks"] > 0:
            log.info(f"Found {n['findpeaks']} findPeaks reports")

        n["tagDir"] = self.homer_tagdirectory()
        if n["tagDir"] > 0:
            log.info(f"Found {n['tagDir']} tagDir reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
