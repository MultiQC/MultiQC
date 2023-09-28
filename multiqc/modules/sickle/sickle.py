""" MultiQC module to parse output from sickle """

import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Sickle",
            anchor="sickle",
            href="https://github.com/najoshi/sickle",
            info="A windowed adaptive trimming tool for FASTQ files using quality.",
            # No DOI // doi=
        )

        # parse list of log files
        self.sickle_data = dict()
        for f in self.find_log_files("sickle"):
            parsed_data = self.parse_logs(f["f"])
            if len(parsed_data):
                if f["s_name"] in self.sickle_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
                self.sickle_data[f["s_name"]] = parsed_data
                self.add_data_source(f)

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, f["s_name"])

        self.sickle_data = self.ignore_samples(self.sickle_data)

        # no file found
        if len(self.sickle_data) == 0:
            raise ModuleNoSamplesFound

        self.write_data_file(self.sickle_data, "multiqc_sickle")

        self.sickle_general_stats_table()
        self.read_count_plot()

    def parse_logs(self, f):
        """Parse the Sickle standard output"""
        regexes = [
            # Paired-end
            ["reads_paired_kept", re.compile("FastQ paired records kept: ([\d,]+) .*")],
            ["reads_single_kept", re.compile("FastQ single records kept: ([\d,]+).*")],
            ["reads_paired_discarded", re.compile("FastQ paired records discarded: ([\d,]+) .*")],
            ["reads_single_discarded", re.compile("FastQ single records discarded: ([\d,]+) .*")],
            # Single-end
            ["reads_single_kept", re.compile("FastQ records kept: ([\d,]+)")],
            ["reads_single_discarded", re.compile("FastQ records discarded: ([\d,]+)")],
        ]
        data = {}
        for l in f.splitlines():
            s = l.split(":")
            # Search regexes for overview stats# Search regexes for overview stats
            for k, regex in regexes:
                match = regex.search(l)
                if match:
                    data[k] = int(match.group(1).replace(",", ""))

        # if paired data : compute total for general stats table
        if "reads_single_kept" in data:
            data["reads_total_kept"] = data.get("reads_paired_kept", 0) + data["reads_single_kept"]

        if "reads_single_discarded" in data:
            data["reads_total_discarded"] = data.get("reads_paired_discarded", 0) + data["reads_single_discarded"]

        # Calculate percentage discarded
        if "reads_total_kept" in data and "reads_single_discarded" in data:
            data["percentage_discarded"] = (
                float(data["reads_total_discarded"]) / float(data["reads_total_kept"] + data["reads_total_discarded"])
            ) * 100.0

        return data

    def sickle_general_stats_table(self):
        """Take the parsed stats from the sickle report and add
        number of kept reads into the generam stats table"""

        headers = {}
        headers["percentage_discarded"] = {
            "title": "% Reads discarded",
            "description": "Percentage of reads discarded",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "OrRd",
        }
        headers["reads_total_kept"] = {
            "title": "{} Reads kept".format(config.read_count_prefix),
            "description": "Number of reads kept ({})".format(config.read_count_desc),
            "modify": lambda x: x * config.read_count_multiplier,
            "scale": "BuGn",
        }
        self.general_stats_addcols(self.sickle_data, headers)

    def read_count_plot(self):
        """Stacked bar plot showing counts of reads"""
        pconfig = {
            "id": "sickle_sequence_counts_plot",
            "title": "Sickle: Sequence Counts",
            "ylab": "Number of reads",
            "cpswitch_counts_label": "Number of reads",
        }
        pcats = OrderedDict()
        pcats["reads_paired_kept"] = {"name": "Paired reads kept", "color": "#1f78b4"}
        pcats["reads_single_kept"] = {"name": "Single reads kept", "color": "#a6cee3"}
        pcats["reads_paired_discarded"] = {"name": "Paired reads discarded", "color": "#ff7f00"}
        pcats["reads_single_discarded"] = {"name": "Single reads discarded", "color": "#fdae61"}

        self.add_section(
            name="Sequence counts",
            anchor="sickle_sequence_counts",
            description="Number of reads discarded by Sickle.",
            plot=bargraph.plot(self.sickle_data, pcats, pconfig),
        )
