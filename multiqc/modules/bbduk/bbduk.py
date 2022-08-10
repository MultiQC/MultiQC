#!/usr/bin/env python
""" Module to parse output from BBDuk """
from __future__ import print_function
from collections import OrderedDict
from multiqc.utils import config
from multiqc.plots import beeswarm
from multiqc.modules.base_module import BaseMultiqcModule

import logging

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """BBDuk Module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="BBDuk",
            anchor="bbduk",
            href="https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/",
            info="""is a tool performing common data-quality-related trimming, 
			filtering, and masking operations with a kmer based approach""",
            ## One publication, but only for the merge tool:
            # doi="10.1371/journal.pone.0185056",
        )

        ## Define the main bbduk multiqc data object
        self.bbduk_data = dict()

        for f in self.find_log_files("bbduk", filehandles=True):
            self.parse_logs(f)

        self.bbduk_data = self.ignore_samples(self.bbduk_data)

        if len(self.bbduk_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.bbduk_data)))

        ## DEBUGGING
        print(self.bbduk_data)

        # Write data to file
        self.write_data_file(self.bbduk_data, "bbduk")

        # self.bbduk_general_stats()
        # self.bbduk_beeswarm_plot()

    def parse_logs(self, logfile):
        """Parses a BBDuk stdout saved in a file"""
        ## Assume take from the name of the file as the includes pair names,
        ## which we can't 'collapse' into a single one.
        s_name = logfile["fn"]
        s_name = self.clean_s_name(s_name, logfile)

        if self.bbduk_data.get(s_name) is not None:
            log.warn("Duplicate sample name found based on filename! Overwriting: {}".format(s_name))

        self.bbduk_data[s_name] = {}
        self.add_data_source(logfile, s_name=s_name)
        file_content = logfile["f"]

        for l in file_content:
            ## Find line after loading reads, and remove suffixes for sample name

            for cat in [
                # "Input",
                "QTrimmed",
                "KTrimmed",
                "Trimmed by overlap",
                "Low quality discards",
                "Low entropy discards",
                "Total Removed",
                "Result",
            ]:
                if cat in l:
                    self.bbduk_data[s_name][cat + " reads"] = int(grab_reads(l))
                    self.bbduk_data[s_name][cat + " percent"] = float(grab_perc(l))


def grab_reads(l):
    """Extracts read counts from STDOUT entry"""
    return l.split(":")[1].lstrip().split(" ")[0]


def grab_perc(l):
    """Extracts percent from STDOUT entry"""
    return l.split(":")[1].lstrip().split(" ")[2].strip("(|)|%")
