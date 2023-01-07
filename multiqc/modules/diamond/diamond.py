""" MultiQC module to parse output from DIAMOND """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="DIAMOND",
            anchor="diamond",
            href="https://github.com/bbuchfink/diamond",
            info="a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data.",
            doi="10.1038/s41592-021-01101-x",
        )

        # Find and load any DIAMOND reports
        self.diamond_data = dict()

        for f in self.find_log_files("diamond", filehandles=True):
            self.parse_logs(f)

        # Filter to strip out ignored sample names
        self.diamond_data = self.ignore_samples(self.diamond_data)

        if len(self.diamond_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.diamond_data)))

        # Write parsed report data to file
        self.write_data_file(self.diamond_data, "diamond")
        self.diamond_general_stats()
        self.diamond_barplot()

    def parse_logs(self, f):
        """Parsing logs""" ""
        for l in f["f"]:
            if "queries aligned" in l:
                self.add_data_source(f)
                if f["s_name"] in self.diamond_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
                self.diamond_data[f["s_name"]] = {"queries_aligned": int(l.split(" ")[0])}

    def diamond_general_stats(self):
        """Diamond General Stats Table"""
        headers = OrderedDict()
        headers["queries_aligned"] = {
            "title": "Queries aligned",
            "description": "number of queries aligned",
            "scale": "YlGn",
        }
        self.general_stats_addcols(self.diamond_data, headers)

    def diamond_barplot(self):
        """Barplot of number of queries aligned"""
        cats = OrderedDict()
        cats["queries_aligned"] = {"name": "Queries Aligned", "color": "#7cb5ec"}
        config = {
            "id": "diamond-barplot",
            "title": "Diamond: Number of queries aligned",
            "ylab": "Number of queries",
        }
        self.add_section(
            name="Queries aligned",
            anchor="barplot",
            description="Shows the number of queries that were aligned to the diamond database.",
            plot=bargraph.plot(self.diamond_data, cats, config),
        )
