import logging
import re

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="telseq",
            anchor="telseq",
            href="https://github.com/zd1/telseq",
            info="Estimates telomere length from whole genome sequencing data (BAMs).",
            extra="""
            Telomeres play a key role in replicative ageing and undergo age-dependent attrition in vivo.
            TelSeq measures average telomere length from whole genome or exome shotgun sequence data.
            """,
            doi="10.1093/nar/gku181",
        )

        # Find and load any telseq reports
        self.telseq_data = dict()

        # Parse the output files
        self.parse_telseq_data()

        # Filter to strip out ignored sample names
        self.telseq_data = self.ignore_samples(self.telseq_data)

        if len(self.telseq_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.telseq_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.telseq_data, "multiqc_telseq")

        # Basic Stats Table
        self.telseq_general_stats_table()


   def parse_telseq_data(self):
        """Go through log file looking for telseq output"""
        for f in self.find_log_files("librarian", filehandles=True):
            headers = f["f"].readline().strip().split("\t")
            for line in f["f"]:
                data = dict(zip(headers, line.strip().split("\t")))
                s_name = self.clean_s_name(data["Sample"], f)
                tel_length = 0.0
                try:
                    tel_length = float(data["LENGTH_ESTIMATE"])
                except Exception:
                    tel_length = data["LENGTH_ESTIMATE"]
                if s_name in self.librarian_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.telseq_data[s_name] = tel_length
                self.add_data_source(f, s_name=s_name)

    def telseq_general_stats_table(self):
        """Take the parsed stats from the telseq report and add it to the
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
        self.general_stats_addcols(self.telseq_data, headers)
