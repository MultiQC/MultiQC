""" MultiQC module to parse output from DIAMOND """

import logging
import re
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

VERSION_REGEX = r"diamond v([\d\.]+)"


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
            raise ModuleNoSamplesFound

        log.info("Found {} reports".format(len(self.diamond_data)))

        # Write parsed report data to file
        self.write_data_file(self.diamond_data, "diamond")
        self.diamond_general_stats()

    def parse_logs(self, f):
        """Parsing logs"""
        s_name = self.clean_s_name(f["root"], f)
        for l in f["f"]:
            # Try to get the sample name from --out, otherwise --query, fallback to directory name
            if "diamond blastx" in l and "--out" in l:
                s_name = l.split("--out ")[1].split(" ")[0]
                s_name = self.clean_s_name(s_name, f)
            elif "diamond blastx" in l and "--query" in l:
                s_name = l.split("--query ")[1].split(" ")[0]
                s_name = self.clean_s_name(s_name, f)

            # Get version
            version_match = re.search(VERSION_REGEX, l)
            if version_match:
                self.add_software_version(version_match.group(1), s_name)

            if "queries aligned" in l:
                self.add_data_source(f)
                if s_name in self.diamond_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.diamond_data[s_name] = {"queries_aligned": int(l.split(" ")[0])}

    def diamond_general_stats(self):
        """Diamond General Stats Table"""
        headers = OrderedDict()
        headers["queries_aligned"] = {
            "title": "Queries aligned",
            "description": "number of queries aligned",
            "scale": "YlGn",
            "format": "{:,.0f}",
        }
        self.general_stats_addcols(self.diamond_data, headers)
