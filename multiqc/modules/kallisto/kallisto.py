""" MultiQC module to parse output from Kallisto """


import logging
import os
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Kallisto module"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Kallisto",
            anchor="kallisto",
            href="http://pachterlab.github.io/kallisto/",
            info="is a program for quantifying abundances of transcripts from RNA-Seq data.",
            doi="10.1038/nbt.3519",
        )

        # Find and load any Kallisto reports
        self.kallisto_data = dict()
        for f in self.find_log_files("kallisto", filehandles=True):
            self.parse_kallisto_log(f)

        # Filter to strip out ignored sample names
        self.kallisto_data = self.ignore_samples(self.kallisto_data)

        if len(self.kallisto_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.kallisto_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.kallisto_data, "multiqc_kallisto")

        # Basic Stats Table
        self.kallisto_general_stats_table()

        # Alignment Rate Plot
        self.add_section(plot=self.kallisto_alignment_plot())

    def parse_kallisto_log(self, f):
        s_name = total_reads = paligned_reads = fraglength = None
        for line in f["f"]:
            # Get input filename
            match = re.search(r"\[quant\] will process (pair|file|sample) 1: (\S+)", line)
            if match:
                s_name = self.clean_s_name(os.path.basename(match.group(2)), f)

            if s_name is not None:
                # Alignment rates
                aligned = re.search(r"\[quant\] processed ([\d,]+) reads, ([\d,]+) reads pseudoaligned", line)
                if aligned:
                    total_reads = float(aligned.group(1).replace(",", ""))
                    paligned_reads = float(aligned.group(2).replace(",", ""))

                # Paired end fragment lengths
                flength = re.search(r"\[quant\] estimated average fragment length: ([\d\.]+)", line)
                if flength:
                    fraglength = float(flength.group(1).replace(",", ""))

                if "quantifying the abundances" in line:
                    if s_name in self.kallisto_data:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name)
                    self.kallisto_data[s_name] = {
                        "total_reads": total_reads,
                        "pseudoaligned_reads": paligned_reads,
                        "not_pseudoaligned_reads": total_reads - paligned_reads,
                    }
                    try:
                        self.kallisto_data[s_name]["percent_aligned"] = (paligned_reads / total_reads) * 100
                    except ZeroDivisionError:
                        self.kallisto_data[s_name]["percent_aligned"] = 0.0
                    if fraglength is not None:
                        self.kallisto_data[s_name]["fragment_length"] = fraglength
                    s_name = total_reads = paligned_reads = fraglength = None

    def kallisto_general_stats_table(self):
        """Take the parsed stats from the Kallisto report and add it to the
        basic stats table at the top of the report"""

        headers = {
            "fragment_length": {
                "title": "Frag Length",
                "description": "Estimated average fragment length",
                "min": 0,
                "suffix": "bp",
                "scale": "RdYlGn",
            },
            "percent_aligned": {
                "title": "% Aligned",
                "description": "% processed reads that were pseudoaligned",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            },
            "pseudoaligned_reads": {
                "title": f"{config.read_count_prefix} Aligned",
                "description": f"Pseudoaligned reads ({config.read_count_desc})",
                "min": 0,
                "scale": "PuRd",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
        }
        self.general_stats_addcols(self.kallisto_data, headers)

    def kallisto_alignment_plot(self):
        # Specify the order of the different possible categories
        keys = {
            "pseudoaligned_reads": {"color": "#437bb1", "name": "Pseudoaligned"},
            "not_pseudoaligned_reads": {"color": "#b1084c", "name": "Not aligned"},
        }

        # Config for the plot
        config = {
            "id": "kallisto_alignment",
            "title": "Kallisto: Alignment Scores",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(self.kallisto_data, keys, config)
