import sys
sys.path.insert(0, '/workspaces/MultiQC/multiqc')
from base_module import BaseMultiqcModule
import re

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="bamutil_clipoverlap",
            anchor="bamutil_clipoverlap",
            href="http://genome.sph.umich.edu/wiki/BamUtil",
            info="bamUtil is a repository that contains several programs that perform operations on SAM/BAM files. This module parses output from the clipOverlap submodule."
        )
        self.bamutil_data = self.parse_logs()
        if len(self.bamutil_data) == 0:
            raise UserWarning("Could not find any bamUtil clipOverlap reports in the specified input files.")
        self.bamutil_general_stats()

    def parse_logs(self):
        data = dict()
        for f in self.find_log_files("bamutil/clipoverlap", filehandles=True):
            s_name = self.clean_s_name(f["s_name"])
            number_overlapping_pairs = None
            avg_ref_bases_overlapped = None
            for line in f["f"]:
                if line.startswith("Number of overlapping pairs:"):
                    number_overlapping_pairs = int(line.strip().split(":")[1])
                elif line.startswith("Average # Reference Bases Overlapped:"):
                    avg_ref_bases_overlapped = float(line.strip().split(":")[1])
            if number_overlapping_pairs is not None and avg_ref_bases_overlapped is not None:
                data[s_name] = {
                    "bamutil_overlapping_pairs": number_overlapping_pairs,
                    "bamutil_avg_ref_bases_overlapped": avg_ref_bases_overlapped
                }
        return data

    def bamutil_general_stats(self):
        headers = {
            "bamutil_overlapping_pairs": {
                "title": "Overlapping Pairs",
                "description": "Number of overlapping pairs (bamUtil clipOverlap)",
                "min": 0,
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "bamutil_avg_ref_bases_overlapped": {
                "title": "Avg Ref Bases Overlapped",
                "description": "Average number of reference bases overlapped (bamUtil clipOverlap)",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Purples",
            },
        }
        self.general_stats_addcols(self.bamutil_data, headers)
