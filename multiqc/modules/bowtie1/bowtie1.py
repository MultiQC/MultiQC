import logging
import re

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Bowtie 1",
            anchor="bowtie1",
            target="Bowtie 1",
            href="http://bowtie-bio.sourceforge.net/",
            info="Ultrafast, memory-efficient short read aligner.",
            doi="10.1186/gb-2009-10-3-r25",
        )

        # Find and load any Bowtie reports
        self.bowtie_data = dict()
        for f in self.find_log_files("bowtie1"):
            self.parse_bowtie_logs(f)

        # Filter to strip out ignored sample names
        self.bowtie_data = self.ignore_samples(self.bowtie_data)

        if len(self.bowtie_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.bowtie_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.bowtie_data, "multiqc_bowtie1")

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bowtie_general_stats_table()

        # Alignment Rate Plot
        self.bowtie_alignment_plot()

    def parse_bowtie_logs(self, f):
        s_name = f["s_name"]
        parsed_data = {}
        regexes = {
            "reads_processed": r"# reads processed:\s+(\d+)",
            "reads_aligned": r"# reads with at least one(?: reported)? alignment:\s+(\d+)",
            "reads_aligned_percentage": r"# reads with at least one(?: reported)? alignment:\s+\d+\s+\(([\d\.]+)%\)",
            "not_aligned": r"# reads that failed to align:\s+(\d+)",
            "not_aligned_percentage": r"# reads that failed to align:\s+\d+\s+\(([\d\.]+)%\)",
            "multimapped": r"# reads with alignments suppressed due to -m:\s+(\d+)",
            "multimapped_percentage": r"# reads with alignments suppressed due to -m:\s+\d+\s+\(([\d\.]+)%\)",
        }

        for line in f["f"].splitlines():
            # Attempt in vain to find original bowtie1 command, logged by another program
            if "bowtie" in line and "q.gz" in line:
                fqmatch = re.search(r"([^\s,]+\.f(ast)?q.gz)", line)
                if fqmatch:
                    s_name = self.clean_s_name(fqmatch.group(1), f)
                    log.debug(f"Found a bowtie command, updating sample name to '{s_name}'")

            # End of log, reset in case there is another in this file
            if "Overall time:" in line:
                if len(parsed_data) > 0:
                    if s_name in self.bowtie_data:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name)
                    self.bowtie_data[s_name] = parsed_data
                s_name = f["s_name"]
                parsed_data = {}

            for k, r in regexes.items():
                match = re.search(r, line)
                if match:
                    parsed_data[k] = float(match.group(1).replace(",", ""))

        if len(parsed_data) > 0:
            if s_name in self.bowtie_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.bowtie_data[s_name] = parsed_data

    def bowtie_general_stats_table(self):
        """Take the parsed stats from the Bowtie report and add it to the
        basic stats table at the top of the report"""

        headers = {
            "reads_aligned_percentage": {
                "title": "% Aligned",
                "description": "% reads with at least one reported alignment",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            },
            "reads_aligned": {
                "title": f"{config.read_count_prefix} Aligned",
                "description": f"reads with at least one reported alignment ({config.read_count_desc})",
                "min": 0,
                "scale": "PuRd",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
        }
        self.general_stats_addcols(self.bowtie_data, headers)

    def bowtie_alignment_plot(self):
        # Specify the order of the different possible categories
        keys = {
            "reads_aligned": {"color": "#8bbc21", "name": "Aligned"},
            "multimapped": {"color": "#2f7ed8", "name": "Multimapped"},
            "not_aligned": {"color": "#0d233a", "name": "Not aligned"},
        }

        # Config for the plot
        config = {
            "id": "bowtie1_alignment",
            "title": "Bowtie 1: Alignment Scores",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        self.add_section(
            description="This plot shows the number of reads aligning to the reference in different ways.",
            helptext="""
            There are 3 possible types of alignment:
            * **Aligned**: Read has only one occurence in the reference genome.
            * **Multimapped**: Read has multiple occurence.
            * **Not aligned**: Read has no occurence.
            """,
            plot=bargraph.plot(self.bowtie_data, keys, config),
        )
