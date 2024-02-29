""" MultiQC submodule to parse output from Rockhopper summary files
https://cs.wellesley.edu/~btjaden/Rockhopper/ """


import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialize the parent object
        super(MultiqcModule, self).__init__(
            name="Rockhopper",
            anchor="rockhopper",
            href="https://cs.wellesley.edu/~btjaden/Rockhopper/",
            info="""
            is a comprehensive and user-friendly system
            for computational analysis of bacterial RNA-seq data.
            It can align reads to genomes, assemble transcripts,
            identify transcript boundaries, and discover novel
            transcripts such as small RNAs.
            """,
            doi=["10.1016/j.ymeth.2019.03.026", "10.1186/s13059-014-0572-2", "10.1093/nar/gkt444"],
        )

        # Set up vars
        self.rh_data = dict()

        # Parse summary file
        for f in self.find_log_files("rockhopper"):
            self.parse_rockhopper_summary(f)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        # Filter to strip out ignored sample names
        self.rh_data = self.ignore_samples(self.rh_data)

        if len(self.rh_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.rh_data)} reports")

        # Write to file
        self.write_data_file(self.rh_data, "multiqc_rockhopper")

        # Basic stats Table
        self.rockhopper_general_stats_table()

        # Alignment bar plot
        self.rockhopper_count_bar_plot()

    def add_results_to_rhdata(self, f, s_name, results):
        """
        Helper function to add parsed results to rhdata
        """
        # Calculate unaligned read count
        total_mapped_reads = sum([v for k, v in results.items() if k != "total-reads"])
        results["unaligned"] = results["total-reads"] - total_mapped_reads

        if s_name in self.rh_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.add_data_source(f, s_name)
        self.rh_data[s_name] = results

    def parse_rockhopper_summary(self, f):
        s_name = None

        # Initialize stats fields
        stats_index = [
            "mRNA-sense",
            "mRNA-antisense",
            "rRNA-sense",
            "rRNA-antisense",
            "tRNA-sense",
            "tRNA-antisense",
            "ncRNA-sense",
            "ncRNA-antisense",
            "unannotated",
        ]

        results = {name: 0 for name in stats_index}
        # Files can have more than one sample in them. Store results for
        # each sample in a dictionary
        results_by_s_name = {}

        # Parse rockhopper output line-by-line since there may be many genomes
        lines = f["f"].split("\n")
        for i, line in enumerate(lines):
            # Get sample name
            if line.startswith("Aligning sequencing reads from file:"):
                # When reaching a new sample add the previous sample to
                # the dictionary of samples along with its results.
                # Then reset the results dictionary for the new sample
                if s_name and s_name not in results_by_s_name:
                    results_by_s_name[s_name] = results
                    results = {name: 0 for name in stats_index}

                s_name = line.split(":", 1)[1].strip()
                s_name = self.clean_s_name(s_name, f)

            elif line.startswith("Aligning sequencing reads from files:"):
                if s_name and s_name not in results_by_s_name:
                    results_by_s_name[s_name] = results
                    results = {name: 0 for name in stats_index}
                s_name = lines[i + 1].strip()
                s_name = self.clean_s_name(s_name, f)

            # Get total number of reads read by rockhopper
            if line.startswith("Total reads:"):
                results["total-reads"] = int(re.search(r"Total reads:\s*(\d*)", line).group(1))

            # Get number of reads aligned to each genome
            elif line.startswith("Successfully aligned reads"):
                # Get number of aligned reads
                genome_reads = int(re.search(r"Successfully aligned reads:\s*(\d*)", line).group(1))

                # Get percent of reads in each category
                stats = [int(re.search(r"(\d+)\%", subline).group(1)) for subline in lines[i + 1 : i + 10]]
                for name, val in zip(stats_index, stats):
                    # Convert percentages to true number of reads in each category
                    results[name] += int(round(val * genome_reads / 100))

        # Make sure the last sample name is added to the results dictionary
        if s_name and s_name not in results_by_s_name:
            results_by_s_name[s_name] = results

        # loop through all samples found and add them to rhdata
        for s_name in results_by_s_name:
            self.add_results_to_rhdata(f, s_name, results_by_s_name[s_name])

    def rockhopper_general_stats_table(self):
        """Take the parsed stats from the Rockhopper summary and add it to the
        basic stats table at the top of the report"""

        headers = {}
        headers["mRNA-sense"] = {
            "title": f"CDS Reads ({config.read_count_prefix})",
            "description": f"Reads aligned to coding regions ({config.read_count_desc})",
            "min": 0,
            "scale": "Blues",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
        }
        headers["mRNA-antisense"] = {
            "title": f"CDS Reads (a/s, {config.read_count_prefix})",
            "description": f"Antisense reads aligned to coding regions ({config.read_count_desc})",
            "min": 0,
            "scale": "Blues",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "hidden": True,
        }
        headers["rRNA-sense"] = {
            "title": f"rRNA Reads ({config.read_count_prefix})",
            "description": f"Reads aligned to rRNA ({config.read_count_desc})",
            "min": 0,
            "scale": "Blues",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
        }
        headers["rRNA-antisense"] = {
            "title": f"rRNA Reads (a/s, {config.read_count_prefix})",
            "description": f"Antisense reads aligned to rRNA ({config.read_count_desc})",
            "min": 0,
            "scale": "Blues",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "hidden": True,
        }
        self.general_stats_addcols(self.rh_data, headers)

    def rockhopper_count_bar_plot(self):
        """Stacked bar plot showing counts of reads"""

        pconfig = {
            "id": "rockhopper_reads_counts_plot",
            "title": "Rockhopper: Alignment types",
            "ylab": "Number of reads",
            "tt_percentage": False,
        }

        # Plot bar graph of groups
        keys = {
            "mRNA-sense": {"name": "mRNA (Sense)"},
            "mRNA-antisense": {"name": "mRNA (Antisense)"},
            "rRNA-sense": {"name": "rRNA (Sense)"},
            "rRNA-antisense": {"name": "rRNA (Antisense)"},
            "tRNA-sense": {"name": "tRNA (Sense)"},
            "tRNA-antisense": {"name": "tRNA (Antisense)"},
            "ncRNA-sense": {"name": "ncRNA (Sense)"},
            "ncRNA-antisense": {"name": "ncRNA (Antisense)"},
            "unannotated": {"name": "Unannotated"},
            "unaligned": {"name": "Unaligned"},
        }

        self.add_section(
            name="Rockhopper",
            anchor="rockhopper_reads_counts",
            description="""
            This plot shows the number of reads mapped to
            different regions of the genome, accounting for
            sense and antisense alignment, if relevant.
            """,
            plot=bargraph.plot(self.rh_data, keys, pconfig),
        )
