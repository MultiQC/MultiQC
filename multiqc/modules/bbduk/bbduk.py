#!/usr/bin/env python
""" Module to parse output from BBDuk """
from __future__ import print_function
from collections import OrderedDict
from multiqc.utils import config
from multiqc.plots import beeswarm, bargraph
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

        # Write data to file
        self.write_data_file(self.bbduk_data, "bbduk")

        self.bbduk_general_stats()
        self.bbduk_bargraph_plot()
        self.bbduk_beeswarm_plot()

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
                elif "Input:" in l:
                    self.bbduk_data[s_name]["Input reads"] = int(grab_reads(l))

    def bbduk_general_stats(self):
        """BBDuk read counts for general stats"""
        headers = OrderedDict()

        headers["Input reads"] = {
            "title": "Total Input Reads ({})".format(config.read_count_prefix),
            "description": "Total number of input reads to BBDuk ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["QTrimmed reads"] = {
            "title": "QTrimmed Reads ({})".format(config.read_count_prefix),
            "description": "QTrimmed Reads removed ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["QTrimmed percent"] = {
            "title": "QTrimmed Percent (%)",
            "description": "Percentage of reads removed through the QTrimmed filter (%)",
            "scale": "Purples",
            "format": "{:,.2f}",
        }
        headers["KTrimmed reads"] = {
            "title": "KTrimmed Reads ({})".format(config.read_count_prefix),
            "description": "KTrimmed Reads removed ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["KTrimmed percent"] = {
            "title": "KTrimmed Percent (%)",
            "description": "Percentage of reads removed through the KTrimmed filter (%)",
            "scale": "Purples",
            "format": "{:,.2f}",
        }
        headers["Trimmed by overlap reads"] = {
            "title": "Trimmed by overlap reads ({})".format(config.read_count_prefix),
            "description": "Reads removed trimmed by overlap ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Trimmed by overlap percent"] = {
            "title": "Trimmed by overlap (%)",
            "description": "Percentage of reads trimmed by overlap (%)",
            "scale": "Purples",
            "format": "{:,.2f}",
        }
        headers["Low quality discards reads"] = {
            "title": "Low quality discards Reads ({})".format(config.read_count_prefix),
            "description": "Reads removed by the Low quality discards filter ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Low quality discards percent"] = {
            "title": "Low quality discards (%)",
            "description": "Percentage of reads removed through the Low quality discards filter (%)",
            "scale": "Purples",
            "format": "{:,.2f}",
        }
        headers["Low entropy discards reads"] = {
            "title": "Low quality discards Reads ({})".format(config.read_count_prefix),
            "description": "Reads removed by the Low entropy discards filter ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Low entropy discards percent"] = {
            "title": "Low entropy discards (%)",
            "description": "Percentage of reads removed through the Low entropy discards filter (%)",
            "scale": "Purples",
            "format": "{:,.2f}",
        }
        headers["Total Removed reads"] = {
            "title": "Total Reads Removed ({})".format(config.read_count_prefix),
            "description": "Total Reads removed ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Total Removed percent"] = {
            "title": "Total Reads Removed (%)",
            "description": "Percentage of reads removed after filtering(%)",
            "scale": "Purples",
            "format": "{:,.2f}",
        }
        headers["Result reads"] = {
            "title": "Remaining Reads ({})".format(config.read_count_prefix),
            "description": "Remaining Reads removed ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Result percent"] = {
            "title": "Total Reads Remaining (%)",
            "description": "Percentage of reads retained after filtering (%)",
            "scale": "Purples",
            "format": "{:,.2f}",
        }

        self.general_stats_addcols(self.bbduk_data, headers)

    def bbduk_bargraph_plot(self):
        """Bargraph summarising reported filtered reads by BBduk"""
        cats = OrderedDict()

        for cat in [
            "Total Removed reads",
            "Result reads",
        ]:
            cats[cat] = {"name": cat}

        config = {
            "id": "bbduk-bargraph",
            "title": "BBDuk filtered reads in percentages",
        }

        self.add_section(
            name="Filtered Reads Percentages",
            anchor="bbduk-bargraph",
            description="Shows summary of reads removed across all BBDuk filters",
            plot=bargraph.plot(self.bbduk_data, cats, config),
        )

    def bbduk_beeswarm_plot(self):
        """
        Beeswarm displaying all possible filtering results reported by BBDuk.

        We don't display this as a barchart as the total across all categories
        of filters reported don't match exactly the total reads remaining (I
        assume there is additional default filtering carried out)
        """
        headers = OrderedDict()

        reads = {
            "min": 0,
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "suffix": "{} reads".format(config.read_count_prefix),
            "decimalPlaces": 0,
            "shared_key": "read_count",
        }
        headers["Input reads"] = dict(reads, title="Input reads")
        headers["QTrimmed reads"] = dict(reads, title="QTtrimmed")
        headers["KTrimmed reads"] = dict(reads, title="KTrimmed")
        headers["Trimmed by overlap reads"] = dict(reads, title="Trimmed by overlap")
        headers["Low quality discards reads"] = dict(reads, title="Low quality discards")
        headers["Low entropy discards reads"] = dict(reads, title="Low entropy discards")
        headers["Total Removed reads"] = dict(reads, title="Total removed")
        headers["Result reads"] = dict(reads, title="Final reads")

        self.add_section(
            name="Filtered Reads",
            anchor="bbduk-beeswarm",
            description="Shows the number of reads removed by various BBDuk filters",
            plot=beeswarm.plot(self.bbduk_data, headers, {"id": "bbduk-beeswarm"}),
        )


## Parsing helper functions


def grab_reads(l):
    """Extracts read counts from STDOUT entry"""
    return l.split(":")[1].lstrip().split(" ")[0]


def grab_perc(l):
    """Extracts percent from STDOUT entry"""
    return l.split(":")[1].lstrip().split(" ")[2].strip("(|)|%")
