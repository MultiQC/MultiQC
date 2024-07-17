import logging
import re

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="leeHom",
            anchor="leehom",
            href="https://github.com/grenaud/leeHom",
            info="Bayesian reconstruction of ancient DNA",
            extra="""
            leeHom is a Bayesian maximum a posteriori algorithm for stripping
            sequencing adapters and merging overlapping portions of reads.
            The algorithm is mostly aimed at ancient DNA and Illumina data but
            can be used for any dataset.
            """,
            doi="10.1093/nar/gku699",
        )

        # Find and load any leeHom reports
        self.leehom_data = dict()

        for f in self.find_log_files("leehom", filehandles=True):
            parsed_data = self.parse_leehom_logs(f)
            if parsed_data is not None and len(parsed_data) > 0:
                self.leehom_data[f["s_name"]] = parsed_data
                self.add_data_source(f, f["s_name"])

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, f["s_name"])

        # Filter to strip out ignored sample names
        self.leehom_data = self.ignore_samples(self.leehom_data)

        if len(self.leehom_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.leehom_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.leehom_data, "multiqc_leehom")

        # Basic Stats Table
        self.leehom_general_stats_table()

    def parse_leehom_logs(self, f):
        """Go through log file looking for leehom output"""
        regexes = {
            "total": r"Total reads[\s\:]+(\d+)",
            "merged_trimming": r"Merged \(trimming\)\s+(\d+)",
            "merged_overlap": r"Merged \(overlap\)\s+(\d+)",
            "kept": r"Kept PE/SR\s+(\d+)",
            "trimmed": r"Trimmed SR\s+(\d+)",
            "adapter_dimers_chimeras": r"Adapter dimers/chimeras\s+(\d+)",
            "failed_key": r"Failed Key\s+(\d+)",
        }
        parsed_data = dict()
        for line in f["f"]:
            # Search regexes for overview stats
            for k, r in regexes.items():
                match = re.search(r, line)
                if match:
                    parsed_data[k] = int(match.group(1))
        return parsed_data

    def leehom_general_stats_table(self):
        """Take the parsed stats from the leeHom report and add it to the
        basic stats table at the top of the report"""

        headers = {}
        headers["merged_trimming"] = {
            "title": f"{config.read_count_prefix} Merged (Trimming)",
            "description": f"Merged clusters from trimming ({config.read_count_desc})",
            "min": 0,
            "scale": "PuRd",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
        }
        headers["merged_overlap"] = {
            "title": f"{config.read_count_prefix} Merged (Overlap)",
            "description": f"Merged clusters from overlapping reads ({config.read_count_desc})",
            "min": 0,
            "scale": "PuRd",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
        }
        self.general_stats_addcols(self.leehom_data, headers)
