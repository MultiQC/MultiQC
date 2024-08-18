"""MultiQC module to parse output from HOMER findpeaks"""

import logging
import os

# Initialise the logger
log = logging.getLogger(__name__)


class FindPeaksReportMixin:
    def parse_homer_findpeaks(self):
        """Find HOMER findpeaks logs and parse their data"""

        self.homer_findpeaks = dict()
        for f in self.find_log_files("homer/findpeaks", filehandles=True):
            self.parse_find_peaks(f)

        # Filter to strip out ignored sample names
        self.homer_findpeaks = self.ignore_samples(self.homer_findpeaks)
        if not self.homer_findpeaks:
            return 0

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.homer_findpeaks, "multiqc_homer_findpeaks")

        # General Stats Table
        stats_headers = {
            "approximate_ip_efficiency": {
                "title": "% Efficiency",
                "description": "Approximate IP efficiency",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "total_peaks": {
                "title": "Total Peaks",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "GnBu",
            },
            "expected_tags_per_peak": {
                "title": "Tags/Peak",
                "description": "Expected tags per peak",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "PuRd",
            },
        }
        self.general_stats_addcols(self.homer_findpeaks, stats_headers, "findpeaks")

        return len(self.homer_findpeaks)

    def parse_find_peaks(self, f):
        """Parse HOMER findPeaks file headers."""
        parsed_data = dict()
        s_name = f["s_name"]
        for line in f["f"]:
            # Start of data
            if line.strip() and not line.strip().startswith("#"):
                break
            # Automatically parse header lines by = symbol
            s = line[2:].split("=")
            if len(s) > 1:
                k = s[0].strip().replace(" ", "_").lower()
                v = s[1].strip().replace("%", "")
                try:
                    parsed_data[k] = float(v)
                except ValueError:
                    parsed_data[k] = v
                if k == "tag_directory":
                    s_name = self.clean_s_name(os.path.basename(v), f, root=os.path.dirname(v))

        if len(parsed_data) > 0:
            if s_name in self.homer_findpeaks:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            self.add_data_source(f, s_name, section="findPeaks")
            self.homer_findpeaks[s_name] = parsed_data
