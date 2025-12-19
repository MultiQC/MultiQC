import logging
import re
from typing import Dict, Optional

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    RiboDetector is a deep learning-based tool for rapid and accurate identification
    of ribosomal RNA sequences from RNA-seq data.

    The module parses log files generated when RiboDetector is run with the `--log` option (v0.3.0+).

    Supported metrics:
    - Total reads processed
    - Predicted rRNA reads
    - Predicted non-rRNA reads
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="RiboDetector",
            anchor="ribodetector",
            href="https://github.com/hzi-bifo/RiboDetector",
            info="Accurate and rapid ribosomal RNA detection based on deep learning.",
            doi="10.1093/nar/gkac112",
        )

        # Find and parse log files
        self.ribodetector: Dict[str, Dict] = {}
        for f in self.find_log_files("ribodetector", filehandles=True):
            parsed = self.parse_ribodetector_log(f)
            if parsed:
                s_name = self.clean_s_name(f["s_name"], f)
                if s_name in self.ribodetector:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.ribodetector[s_name] = parsed
                self.add_data_source(f, s_name)

        # Filter ignored samples
        self.ribodetector = self.ignore_samples(self.ribodetector)

        if len(self.ribodetector) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.ribodetector)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add to General Stats table
        self.ribodetector_general_stats()

        # Add rRNA detection plot section
        self.add_section(
            name="rRNA Detection",
            anchor="ribodetector_rrna",
            description="Classification of reads as rRNA or non-rRNA by RiboDetector.",
            plot=self.ribodetector_bargraph(),
        )

        # Write parsed data to file
        self.write_data_file(self.ribodetector, "multiqc_ribodetector")

    def parse_ribodetector_log(self, f) -> Optional[Dict]:
        """Parse a RiboDetector log file."""
        data: Dict = {}

        # Regex to strip ANSI escape codes
        ansi_escape = re.compile(r"\x1b\[[0-9;]*m")

        for line in f["f"]:
            # Strip ANSI color codes
            clean_line = ansi_escape.sub("", line)

            # Parse total sequences processed
            # Format: "Processed 49593 sequences in total"
            m = re.search(r"Processed\s+([\d,]+)\s+sequences in total", clean_line, re.IGNORECASE)
            if m:
                data["total_reads"] = int(m.group(1).replace(",", ""))
                continue

            # Parse non-rRNA sequences
            # Format: "Detected 49587 non-rRNA sequences"
            m = re.search(r"Detected\s+([\d,]+)\s+non-rRNA sequences", clean_line, re.IGNORECASE)
            if m:
                data["nonrrna_reads"] = int(m.group(1).replace(",", ""))
                continue

            # Parse rRNA sequences
            # Format: "Detected 6 rRNA sequences"
            m = re.search(r"Detected\s+([\d,]+)\s+rRNA sequences", clean_line, re.IGNORECASE)
            if m:
                data["rrna_reads"] = int(m.group(1).replace(",", ""))
                continue

        # Calculate percentages if we have the data
        if "total_reads" in data and data["total_reads"] > 0:
            if "rrna_reads" in data:
                data["percent_rrna"] = (data["rrna_reads"] / data["total_reads"]) * 100
            if "nonrrna_reads" in data:
                data["percent_nonrrna"] = (data["nonrrna_reads"] / data["total_reads"]) * 100

        return data if data else None

    def ribodetector_general_stats(self):
        """Add columns to the General Statistics table."""
        headers = {
            "percent_rrna": {
                "title": "% rRNA",
                "description": "Percentage of reads classified as rRNA",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "OrRd",
            },
            "nonrrna_reads": {
                "title": f"{config.read_count_prefix} Non-rRNA",
                "description": f"Reads classified as non-rRNA ({config.read_count_desc})",
                "min": 0,
                "scale": "Greens",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "hidden": True,
            },
            "rrna_reads": {
                "title": f"{config.read_count_prefix} rRNA",
                "description": f"Reads classified as rRNA ({config.read_count_desc})",
                "min": 0,
                "scale": "Reds",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.ribodetector, headers)

    def ribodetector_bargraph(self):
        """Create a bar graph showing rRNA vs non-rRNA reads."""
        keys = {
            "nonrrna_reads": {"color": "#437bb1", "name": "Non-rRNA"},
            "rrna_reads": {"color": "#e63946", "name": "rRNA"},
        }

        pconfig = {
            "id": "ribodetector_classification",
            "title": "RiboDetector: Read Classification",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(self.ribodetector, keys, pconfig)
