""" MultiQC module to parse output from Tophat """


import logging
import os
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Tophat",
            anchor="tophat",
            href="https://ccb.jhu.edu/software/tophat/",
            info="is a fast splice junction mapper for RNA-Seq reads. "
            "It aligns RNA-Seq reads to mammalian-sized genomes.",
            doi=["10.1186/gb-2013-14-4-r36", "10.1093/bioinformatics/btp120"],
        )

        # Find and load any Tophat reports
        self.tophat_data = dict()
        for f in self.find_log_files("tophat"):
            parsed_data = self.parse_tophat_log(f["f"])
            if parsed_data is not None:
                if f["s_name"] == "align" or f["s_name"] == "align_summary.txt":
                    s_name = os.path.basename(f["root"])
                else:
                    s_name = f["s_name"].split("align_summary.txt", 1)[0]
                s_name = self.clean_s_name(s_name, f)
                if s_name in self.tophat_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name)
                self.tophat_data[s_name] = parsed_data

        # Filter to strip out ignored sample names
        self.tophat_data = self.ignore_samples(self.tophat_data)

        if len(self.tophat_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.tophat_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.tophat_data, "multiqc_tophat.txt")

        # Basic Stats Table
        self.tophat_general_stats_table()

        # Alignment Rate Plot
        self.tophat_alignment_plot()

    def parse_tophat_log(self, raw_data):
        """Parse the Tophat alignment log file."""

        if "Aligned pairs" in raw_data:
            # Paired end data
            regexes = {
                "overall_aligned_percent": r"([\d\.]+)% overall read mapping rate.",
                "concordant_aligned_percent": r"([\d\.]+)% concordant pair alignment rate.",
                "aligned_total": r"Aligned pairs:\s+(\d+)",
                "aligned_multimap": r"Aligned pairs:\s+\d+\n\s+of these:\s+(\d+)",
                "aligned_discordant": r"(\d+) \([\s\d\.]+%\) are discordant alignments",
                "total_reads": r"[Rr]eads:\n\s+Input\s*:\s+(\d+)",
            }
        else:
            # Single end data
            regexes = {
                "total_reads": r"[Rr]eads:\n\s+Input\s*:\s+(\d+)",
                "aligned_total": r"Mapped\s*:\s+(\d+)",
                "aligned_multimap": r"of these\s*:\s+(\d+)",
                "overall_aligned_percent": r"([\d\.]+)% overall read mapping rate.",
            }

        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = float(r_search.group(1))

        # Exit if we didn't manage to parse enough fields - probably not a TopHat log
        # Note that Bowtie2 / HiSAT2 logs contain some but not all of these strings
        if len(parsed_data) < 4:
            return None

        parsed_data["concordant_aligned_percent"] = parsed_data.get("concordant_aligned_percent", 0)
        parsed_data["aligned_total"] = parsed_data.get("aligned_total", 0)
        parsed_data["aligned_multimap"] = parsed_data.get("aligned_multimap", 0)
        parsed_data["aligned_discordant"] = parsed_data.get("aligned_discordant", 0)
        parsed_data["unaligned_total"] = parsed_data["total_reads"] - parsed_data["aligned_total"]
        parsed_data["aligned_not_multimapped_discordant"] = (
            parsed_data["aligned_total"] - parsed_data["aligned_multimap"] - parsed_data["aligned_discordant"]
        )
        return parsed_data

    def tophat_general_stats_table(self):
        """Take the parsed stats from the Tophat report and add it to the
        basic stats table at the top of the report"""

        headers = {
            "overall_aligned_percent": {
                "title": "% Aligned",
                "description": "overall read mapping rate",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            },
            "aligned_not_multimapped_discordant": {
                "title": f"{config.read_count_prefix} Aligned",
                "description": f"Aligned reads, not multimapped or discordant ({config.read_count_desc})",
                "min": 0,
                "scale": "PuRd",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
        }
        self.general_stats_addcols(self.tophat_data, headers)

    def tophat_alignment_plot(self):
        # Specify the order of the different possible categories
        keys = {
            "aligned_not_multimapped_discordant": {"color": "#437bb1", "name": "Aligned"},
            "aligned_multimap": {"color": "#f7a35c", "name": "Multimapped"},
            "aligned_discordant": {"color": "#e63491", "name": "Discordant mappings"},
            "unaligned_total": {"color": "#7f0000", "name": "Not aligned"},
        }

        # Config for the plot
        config = {
            "id": "tophat_alignment",
            "title": "Tophat: Alignment Scores",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        self.add_section(plot=bargraph.plot(self.tophat_data, keys, config))
