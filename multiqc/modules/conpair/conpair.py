""" MultiQC module to parse output from Conpair """


import logging
import re

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Conpair module class.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Conpair",
            anchor="conpair",
            href="https://github.com/nygenome/Conpair",
            info="is a fast and robust method dedicated for human tumor-normal "
            "studies to perform concordance verification, as well as "
            "cross-individual contamination level estimation in "
            "whole-genome and whole-exome sequencing experiments.",
            doi="bioinformatics/btw389",
        )

        self.conpair_data = dict()

        for f in self.find_log_files("conpair/concordance"):
            self.parse_conpair_logs(f)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        for f in self.find_log_files("conpair/contamination"):
            self.parse_conpair_logs(f)

        # Filter to strip out ignored sample names
        self.conpair_data = self.ignore_samples(self.conpair_data)

        if len(self.conpair_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.conpair_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.conpair_data, "multiqc_conpair")

        # Basic Stats Table
        self.conpair_general_stats_table()

    def parse_conpair_logs(self, f):
        """Go through log file looking for conpair concordance or contamination output
        One parser to rule them all."""

        conpair_regexes = {
            "concordance_concordance": [r"Concordance: ([\d\.]+)%\nBased.*", r"(\d\.\d+)\nBased.*"],
            "concordance_used_markers": [r"Based on (\d+)/\d+ markers"],
            "concordance_total_markers": [r"Based on \d+/(\d+) markers"],
            "concordance_marker_threshold": [r"\(coverage per marker threshold : (\d+) reads\)"],
            "concordance_min_mapping_quality": [r"Minimum mappinq quality: (\d+)"],
            "concordance_min_base_quality": [r"Minimum base quality: (\d+)"],
            "contamination_normal": [r"Normal sample contamination level: ([\d\.]+)%"],
            "contamination_tumor": [r"Tumor sample contamination level: ([\d\.]+)%"],
        }

        parsed_data = {}
        for k, r_arr in conpair_regexes.items():
            for r in r_arr:
                match = re.search(r, f["f"])
                if match:
                    parsed_data[k] = float(match.group(1))
                    if k == "concordance_concordance" and "Concordance" not in r:
                        parsed_data[k] = 100.0 * float(parsed_data[k])
                    break

        def _cp_type(data):
            if "concordance_concordance" in parsed_data:
                return "concordance"
            elif "contamination_normal" in parsed_data:
                return "contamination"

        if len(parsed_data) > 0:
            if f["s_name"] in self.conpair_data:
                if _cp_type(self.conpair_data[f["s_name"]]) == _cp_type(parsed_data):
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            else:
                self.conpair_data[f["s_name"]] = dict()
            self.add_data_source(f, section=_cp_type(parsed_data))
            self.conpair_data[f["s_name"]].update(parsed_data)

    def conpair_general_stats_table(self):
        """Take the parsed stats from the Conpair report and add it to the
        basic stats table at the top of the report"""

        headers = {}
        headers["concordance_concordance"] = {
            "title": "Concordance",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "RdYlGn",
        }
        headers["contamination_normal"] = {
            "title": "N Contamination",
            "description": "Normal sample contamination level",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "format": "{:,.3f}",
            "scale": "RdYlBu-rev",
        }
        headers["contamination_tumor"] = {
            "title": "T Contamination",
            "description": "Tumor sample contamination level",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "format": "{:,.3f}",
            "scale": "RdYlBu-rev",
        }
        self.general_stats_addcols(self.conpair_data, headers)
