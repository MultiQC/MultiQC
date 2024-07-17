import logging
import re
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc import config

log = logging.getLogger(__name__)

VERSION_REGEX = r"Version ([\d\.]+)"


class MultiqcModule(BaseMultiqcModule):
    """
    The module produces summary statistics from the stdout logging information from the BBDuk tool of the
    [BBTools](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/) suite of tools.

    "Duk" stands for Decontamination Using Kmers. BBDuk was developed to combine
    most common data-quality-related trimming, filtering, and masking operations
    into a single high-performance tool.

    The module can summarise data from the following BBDuk funtionality
    (descriptions from command line help output):

    - `entropy` - entropy filtering
    - `ktrim` - kmer trimming
    - `qtrim` - quality trimming
    - `maq` - read quality filtering
    - `ref` contaminant filtering

    Additional information on the BBMap tools is available on
    [SeqAnswers](http://seqanswers.com/forums/showthread.php?t=41057).
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="BBDuk",
            anchor="bbduk",
            href="https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/",
            info="Common data-quality-related trimming, filtering, and masking operations with a kmer based approach",
            # One publication, but only for the merge tool:
            # doi="10.1371/journal.pone.0185056",
        )

        # Define the main bbduk multiqc data object
        self.bbduk_data: Dict = dict()

        for f in self.find_log_files("bbduk", filehandles=True):
            self.parse_logs(f)

        self.bbduk_data = self.ignore_samples(self.bbduk_data)

        if len(self.bbduk_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.bbduk_data)} reports")

        # Write data to file
        self.write_data_file(self.bbduk_data, "bbduk")

        self.bbduk_general_stats()
        self.bbduk_bargraph_plot()

    def parse_logs(self, f):
        """Parses a BBDuk stdout saved in a file"""

        s_name = f["s_name"]
        for line in f["f"]:
            if "jgi.BBDuk" in line and "in1=" in line:
                s_name = line.split("in1=")[1].split(" ")[0]
                s_name = self.clean_s_name(s_name, f)

            if line.startswith("Version"):
                version_match = re.search(VERSION_REGEX, line)
                if version_match:
                    self.add_software_version(version_match.group(1), s_name)

            if "Input:" in line:
                matches = re.search(r"Input:\s+(\d+) reads\s+(\d+) bases", line)
                if matches:
                    self.add_data_source(f, s_name)
                    if s_name in self.bbduk_data:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.bbduk_data[s_name] = dict()

                    self.bbduk_data[s_name]["Input reads"] = int(matches.group(1))
                    self.bbduk_data[s_name]["Input bases"] = int(matches.group(2))
            # Don't start using regexes until we're in that block
            elif "Input reads" in self.bbduk_data.get(s_name, {}):
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
                    matches = re.search(rf"{cat}:\s+(\d+) reads \(([\d\.]+)%\)\s+(\d+) bases \(([\d\.]+)%\)", line)
                    if matches:
                        self.bbduk_data[s_name][cat + " reads"] = int(matches.group(1))
                        self.bbduk_data[s_name][cat + " reads percent"] = float(matches.group(2))
                        self.bbduk_data[s_name][cat + " bases"] = int(matches.group(3))
                        self.bbduk_data[s_name][cat + " bases percent"] = float(matches.group(4))
                        break
            elif "Reads Processed:" in line:
                return

    def bbduk_general_stats(self):
        """BBDuk read counts for general stats"""
        headers = {
            "Total Removed bases percent": {
                "title": "Bases Removed (%)",
                "description": "Percentage of bases removed after filtering",
                "scale": "YlOrBr",
                "max": 100,
            },
            "Total Removed bases": {
                "title": f"Bases Removed ({config.base_count_prefix})",
                "description": f"Total Bases removed ({config.base_count_desc})",
                "scale": "Reds",
                "shared_key": "base_count",
                "modify": lambda x: x * config.base_count_multiplier,
                "hidden": True,
            },
            "Total Removed reads percent": {
                "title": "Reads Removed (%)",
                "description": "Percentage of reads removed after filtering",
                "scale": "OrRd",
                "max": 100,
            },
            "Total Removed reads": {
                "title": f"Reads Removed ({config.read_count_prefix})",
                "description": f"Total Reads removed ({config.read_count_desc})",
                "scale": "Reds",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
            "Input reads": {
                "title": f"Total Input Reads ({config.read_count_prefix})",
                "description": f"Total number of input reads to BBDuk ({config.read_count_desc})",
                "scale": "Greens",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
        }

        self.general_stats_addcols(self.bbduk_data, headers)

    def bbduk_bargraph_plot(self):
        """
        Violin displaying all possible filtering results reported by BBDuk.

        We don't display this as a barchart as the total across all categories
        of filters reported don't match exactly the total reads remaining (I
        assume there is additional default filtering carried out)
        """
        cats = [
            "Result",
            "QTrimmed",
            "KTrimmed",
            "Trimmed by overlap",
            "Low quality discards",
            "Low entropy discards",
        ]
        pconfig = {
            "id": "bbduk-filtered-barplot",
            "title": "BBDuk: Filtered reads",
            "ylab": "Number of Reads",
            "data_labels": [
                {"name": "Reads", "ylab": "Number of Reads"},
                {"name": "Bases", "ylab": "Number of Base Pairs"},
            ],
        }

        self.add_section(
            name="BBDuk: Filtered Reads",
            anchor="bbduk-filtered-reads",
            description="The number of reads removed by various BBDuk filters",
            plot=bargraph.plot(
                [self.bbduk_data, self.bbduk_data],
                [
                    [f"{cat} reads" for cat in cats],
                    [f"{cat} bases" for cat in cats],
                ],
                pconfig,
            ),
        )
