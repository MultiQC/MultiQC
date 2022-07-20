from multiqc.modules.base_module import BaseMultiqcModule
import logging
import re
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Porechop",
            anchor="porechop",
            href="https://github.com/rrwick/Porechop",
            info="A tool for finding and removing adapters from Oxford Nanopore reads.",
        )

        # Find and load reports
        self.porechop_data = dict()

        # Find all files for porechop
        for f in self.find_log_files("porechop", filehandles=True):
            self.parse_logs(f)

        self.porechop_data = self.ignore_samples(self.porechop_data)

        if len(self.porechop_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.porechop_data)))

        # Write data to file
        self.write_data_file(self.porechop_data, "porechop")

        ## DEBUG
        print(self.porechop_data)

        self.porechop_general_stats()

    def parse_logs(self, logfile):
        """Parsing Logs. Note: careful of ANSI formatting log"""
        file_content = logfile["f"]
        for l in file_content:
            ## Find line after loading reads, and remove suffixes for sample name
            if "Loading reads" in l:
                s_name = next(file_content).rstrip()
                s_name = self.clean_s_name(s_name, logfile)
                self.add_data_source(logfile, s_name=s_name)
                self.porechop_data[s_name] = {}
            ## Find each valid metric, clean up for plain integer
            elif "reads loaded" in l:
                self.porechop_data[s_name]["Input Reads"] = {}
                self.porechop_data[s_name]["Input Reads"] = float(l.split(" ")[0])
            elif "from their start" in l:
                self.porechop_data[s_name]["Start Trimmed"] = {}
                self.porechop_data[s_name]["Start Trimmed"] = float(l.split(" ")[0])
                self.porechop_data[s_name]["Start Trimmed Total"] = float(l.split(" ")[2])
                self.porechop_data[s_name]["Start Trimmed (bp)"] = float(l.split(" ")[10].strip("(").replace(",", ""))
                self.porechop_data[s_name]["Start Trimmed (%)"] = (
                    self.porechop_data[s_name]["Start Trimmed"]
                    / self.porechop_data[s_name]["Start Trimmed Total"]
                    * 100
                )
            elif "from their end" in l:
                self.porechop_data[s_name]["End Trimmed"] = {}
                self.porechop_data[s_name]["End Trimmed"] = float(l.split(" ")[0])
                self.porechop_data[s_name]["End Trimmed Total"] = float(l.split(" ")[2])
                self.porechop_data[s_name]["End Trimmed (bp)"] = float(l.split(" ")[10].strip("(").replace(",", ""))
                self.porechop_data[s_name]["End Trimmed (%)"] = (
                    self.porechop_data[s_name]["End Trimmed"] / self.porechop_data[s_name]["Start Trimmed Total"] * 100
                )
            elif "split based on" in l:
                self.porechop_data[s_name]["Middle Split"] = {}
                self.porechop_data[s_name]["Middle Split"] = float(l.split(" ")[0])
                self.porechop_data[s_name]["Middle Split Total"] = float(l.split(" ")[2])
                self.porechop_data[s_name]["Middle Split (%)"] = (
                    self.porechop_data[s_name]["Middle Split"] / self.porechop_data[s_name]["Middle Split Total"] * 100
                )

    def porechop_general_stats(self):
        """Porechop General Stats Table"""
        headers = OrderedDict()
        headers["Input Reads"] = {
            "title": "Input Reads ({})".format(config.read_count_prefix),
            "description": "Number of reads loaded into Porechop ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["Start Trimmed"] = {
            "title": "Start Trimmed ({})".format(config.read_count_prefix),
            "description": "Number of reads that had adapters trimmed from the start ({})".format(
                config.read_count_prefix
            ),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Start Trimmed (%)"] = {
            "title": "Start Trimmed",
            "description": "Percent of reads that had adapters trimmed from the start",
            "suffix": "%",
            "max": 100,
            "scale": "RdYlGn",
        }
        headers["End Trimmed"] = {
            "title": "End Trimmed ({})".format(config.read_count_prefix),
            "description": "Number of reads that had adapters trimmed from the end ({})".format(
                config.read_count_prefix
            ),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["End Trimmed (%)"] = {
            "title": "End Trimmed",
            "description": "Percent of reads that had adapters trimmed from the end",
            "suffix": "%",
            "max": 100,
            "scale": "RdYlGn",
        }
        headers["Middle Split"] = {
            "title": "Middle Split ({})".format(config.read_count_prefix),
            "description": "Number of reads split based on middle adapters ({})".format(config.read_count_prefix),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Middle Split (%)"] = {
            "title": "Middle Split",
            "description": "Percent of reads that were split based on middle adapters",
            "suffix": "%",
            "max": 100,
            "scale": "RdYlGn",
        }
        # headers["Start Trimmed (%)"]
        # headers["End Trimmed (%)"]
        # headers["Middle Split"]
        self.general_stats_addcols(self.porechop_data, headers)
