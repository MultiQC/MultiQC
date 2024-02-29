""" MultiQC module to parse output from MALT """


import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

VERSION_REGEX = r"Version\s+MALT \(version ([\d\.]+),.*"


class MultiqcModule(BaseMultiqcModule):
    """Malt Module"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="MALT",
            anchor="malt",
            href="http://ab.inf.uni-tuebingen.de/software/malt/",
            info="performs alignment of metagenomic reads against a database of reference sequences (such as NR, GenBank or Silva) and produces a MEGAN RMA file as output.",
            doi="10.1101/050559 ",
        )

        # Find and load Malt reports
        self.malt_data = dict()

        for f in self.find_log_files("malt", filehandles=True):
            self.parse_logs(f)

        self.malt_data = self.ignore_samples(self.malt_data)

        if len(self.malt_data) == 0:
            raise ModuleNoSamplesFound

        # Write data to file
        self.write_data_file(self.malt_data, "malt")

        log.info(f"Found {len(self.malt_data)} reports")

        self.malt_general_stats()
        self.mappability_barplot()
        self.taxonomic_assignation_barplot()

    def parse_logs(self, f):
        """Parses a Malt log file"""
        reading = False
        keys = [
            "Total reads",
            "Assig. Taxonomy",
            "Num. of queries",
            "Aligned queries",
            "Num. alignments",
        ]
        version = None
        for line in f["f"]:
            line = line.rstrip()
            if line.startswith("Version"):
                version_match = re.search(VERSION_REGEX, line)
                if version_match:
                    version = version_match.group(1)

            if line.startswith("+++++ Aligning file:") and reading is False:
                reading = True
                s_name = line.split()[-1]
                s_name = self.clean_s_name(s_name, f)
                if s_name in self.malt_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name=s_name)
                self.malt_data[s_name] = {}
                if version is not None:
                    self.add_software_version(version, s_name)

            elif reading:
                for k in keys:
                    if line.startswith(k):
                        self.malt_data[s_name][k] = int(line.split()[-1].replace(",", ""))
                        if k == "Num. alignments":
                            try:
                                self.malt_data[s_name]["Non mapped"] = (
                                    self.malt_data[s_name]["Num. of queries"] - self.malt_data[s_name]["Total reads"]
                                )
                                self.malt_data[s_name]["No Assig. Taxonomy"] = (
                                    self.malt_data[s_name]["Total reads"] - self.malt_data[s_name]["Assig. Taxonomy"]
                                )
                                try:
                                    self.malt_data[s_name]["Mappability"] = (
                                        float(self.malt_data[s_name]["Total reads"])
                                        / float(self.malt_data[s_name]["Num. of queries"])
                                    ) * 100.0
                                except ZeroDivisionError:
                                    self.malt_data[s_name]["Mappability"] = 0
                                try:
                                    self.malt_data[s_name]["Taxonomic assignment success"] = (
                                        float(self.malt_data[s_name]["Assig. Taxonomy"])
                                        / float(self.malt_data[s_name]["Total reads"])
                                    ) * 100.0
                                except ZeroDivisionError:
                                    self.malt_data[s_name]["Taxonomic assignment success"] = 0
                            except KeyError:
                                pass
                            reading = False

    def mappability_barplot(self):
        """Mappability barplot

        Mappability = (Total reads / Num. of queries) * 100
        """
        cats = {
            "Total reads": {"name": "Mapped reads"},
            "Non mapped": {"name": "Non Mapped reads"},
        }
        config = {
            "id": "malt-mappability-plot",
            "title": "MALT: Metagenomic Mappability",
            "ylab": "Read Counts",
        }
        self.add_section(
            name="Metagenomic Mappability",
            anchor="malt-mappability",
            description="Number of mapped reads.",
            plot=bargraph.plot(self.malt_data, cats, config),
        )

    def taxonomic_assignation_barplot(self):
        """Taxonomic assignment barplot

        Taxonomic assignment success = (Assig. Taxonomy / Total reads) * 100
        """
        cats = ["Assig. Taxonomy", "No Assig. Taxonomy"]
        config = {
            "id": "malt-taxonomic-success-plot",
            "title": "MALT: Taxonomic assignment success",
            "ylab": "Read Counts",
        }
        self.add_section(
            name="Taxonomic assignment success",
            anchor="malt-taxonomic-success",
            description="Shows the number of mapped reads assigned to a taxonomic node.",
            plot=bargraph.plot(self.malt_data, cats, config),
        )

    def malt_general_stats(self):
        """MALT General Statistics table"""
        headers = {
            "Taxonomic assignment success": {
                "title": "% Tax assigned",
                "description": "Percentage of mapped reads assigned to a taxonomic node",
                "suffix": "%",
                "max": 100,
                "scale": "RdYlGn",
            },
            "Assig. Taxonomy": {
                "title": f"{config.read_count_prefix} Tax assigned",
                "description": f"Number of reads assigned to a Taxonomic node ({config.read_count_desc})",
                "scale": "Greens",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
            "Mappability": {
                "title": "% Metagenomic Mapped",
                "description": "Percentage of mapped reads",
                "suffix": "%",
                "max": 100,
                "scale": "RdYlGn",
            },
            "Total reads": {
                "title": f"{config.read_count_prefix} Mapped",
                "description": f"Number of mapped reads ({config.read_count_desc})",
                "scale": "PuBu",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "Num. of queries": {
                "title": f"{config.read_count_prefix} Reads",
                "description": f"Number of reads in sample ({config.read_count_desc})",
                "scale": "Purples",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
        }
        self.general_stats_addcols(self.malt_data, headers)
