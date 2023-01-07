""" Module to parse output from BBDuk """

import logging
import re
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph
from multiqc.utils import config

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

    def parse_logs(self, f):
        """Parses a BBDuk stdout saved in a file"""

        ## Assume take from the name of the file as the includes pair names,
        ## which we can't 'collapse' into a single one.
        if f["s_name"] in self.bbduk_data:
            log.warn("Duplicate sample name found based on filename! Overwriting: {}".format(f["s_name"]))
        self.bbduk_data[f["s_name"]] = {}

        self.add_data_source(f)

        for l in f["f"]:
            ## Find line after loading reads, and remove suffixes for sample name

            if "Input:" in l:
                matches = re.search(r"Input:\s+(\d+) reads\s+(\d+) bases", l)
                if matches:
                    self.bbduk_data[f["s_name"]]["Input reads"] = int(matches.group(1))
                    self.bbduk_data[f["s_name"]]["Input bases"] = int(matches.group(2))
            # Don't start using regexes until we're in that block
            elif "Input reads" in self.bbduk_data[f["s_name"]]:
                cats = [
                    "QTrimmed",
                    "KTrimmed",
                    "Trimmed by overlap",
                    "Low quality discards",
                    "Low entropy discards",
                    "Total Removed",
                    "Result",
                ]
                for cat in cats:
                    matches = re.search(f"{cat}:\s+(\d+) reads \(([\d\.]+)%\)\s+(\d+) bases \(([\d\.]+)%\)", l)
                    if matches:
                        self.bbduk_data[f["s_name"]][cat + " reads"] = int(matches.group(1))
                        self.bbduk_data[f["s_name"]][cat + " percent"] = float(matches.group(2))
                        self.bbduk_data[f["s_name"]][cat + " bases"] = int(matches.group(3))
                        self.bbduk_data[f["s_name"]][cat + " bases percent"] = float(matches.group(4))
                        break
            elif "Reads Processed:" in l:
                return

    def bbduk_general_stats(self):
        """BBDuk read counts for general stats"""
        headers = OrderedDict()

        headers["Total Removed percent"] = {
            "title": "Total Reads Removed (%)",
            "description": "Percentage of reads removed after filtering(%)",
            "scale": "OrRd",
            "max": 100,
        }
        headers["Total Removed reads"] = {
            "title": "Total Reads Removed ({})".format(config.read_count_prefix),
            "description": "Total Reads removed ({})".format(config.read_count_desc),
            "scale": "Reds",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Input reads"] = {
            "title": "Total Input Reads ({})".format(config.read_count_prefix),
            "description": "Total number of input reads to BBDuk ({})".format(config.read_count_desc),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        self.general_stats_addcols(self.bbduk_data, headers)

    def bbduk_bargraph_plot(self):
        """
        Beeswarm displaying all possible filtering results reported by BBDuk.

        We don't display this as a barchart as the total across all categories
        of filters reported don't match exactly the total reads remaining (I
        assume there is additional default filtering carried out)
        """
        cats = [
            "Result reads",
            "QTrimmed reads",
            "KTrimmed reads",
            "Trimmed by overlap reads",
            "Low quality discards reads",
            "Low entropy discards reads",
        ]
        pconfig = {
            "id": "bbduk-filtered-barplot",
            "title": "BBDuk: Percentage Summary",
            "ylab": "Percentage of Reads",
        }

        self.add_section(
            name="BBDuk: Filtered Reads",
            anchor="bbduk-filtered-reads",
            description="The number of reads removed by various BBDuk filters",
            plot=bargraph.plot(
                self.bbduk_data,
                cats,
                pconfig,
            ),
        )
