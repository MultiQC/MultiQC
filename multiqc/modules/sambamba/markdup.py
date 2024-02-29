import logging
import re

from multiqc.plots import bargraph

log = logging.getLogger(__name__)

VERSION_REGEX = r"sambamba ([\d\.]+)"


class SambambaMarkdupMixin:

    """Find and parse Sambamba Markdup output log files and calculate duplication rate"""

    def parse_sambamba_markdup(self):
        # Clean sample name from 'markdup_sample_1' to 'sample_1'
        # Find and load sambamba logs to markdup_data.
        # Regex for key phrases and calculate duplicate rate.
        # Detect if sample is single or paired-end reads. Process differently.
        # Give user warning if redundant samples are found.

        self.markdup_data = dict()

        for f in self.find_log_files("sambamba/markdup"):
            if f["s_name"] in self.markdup_data:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {f['s_name']}")

            # parse sambamba output by sample name
            parsed = self.parse_markdup_stats(f)
            if parsed is not None:
                self.markdup_data[f["s_name"]] = parsed

            # filter away samples if MultiQC user does not want them
            self.markdup_data = self.ignore_samples(self.markdup_data)

            # add sambamba version to software table
            version_match = re.search(VERSION_REGEX, f["f"])
            if version_match:
                self.add_software_version(version_match.group(1), f["s_name"])

            # add results to multiqc_sources.txt
            self.add_data_source(f)

            # warn user if duplicate samples found
            if f["s_name"] in self.markdup_data:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")

        if len(self.markdup_data) == 0:
            return 0

        # Write parsed files to a file
        self.write_data_file(self.markdup_data, "multiqc_markdup", sort_cols=True)

        # Add sambamba markdup to general statistics table
        self.markdup_general_stats_table()

        # Add sambamba markdup bargraph section
        self.markdup_section()

        # return the length of log files parsed
        return len(self.markdup_data)

    def parse_markdup_stats(self, f):
        """
        Input actual content of sambamba markdup log output.
        Regex for important phrases and extract reads.
        Detect and calculate duplicate rate by single/paired-end samples.
        Outputs dictionary of markdup stats of a sample by single/paired-end type.
        This dict is associated with corresponding sample in key of markdup_data.

        NOTE: reads seem to sum to between 98 - 99% of input BAM file. Is Sambamba Markdup not recognizing the entire file?
        """

        regexes = {
            "sorted_end_pairs": r"sorted (\d+) end pairs",
            "single_ends": r"and (\d+) single ends",
            "single_unmatched_pairs": r"among them (\d+) unmatched",
            "duplicate_reads": r"found (\d+) duplicates",
        }
        d = {}
        for key, regex in regexes.items():
            m = re.search(regex, f["f"])
            if m:
                d[key] = int(m[1])
        if len(d) != len(regexes):
            log.debug(f"Could not find all markdup fields for '{f['fn']}' - skipping")
            return None

        # Calculate duplicate rate
        # NB: Single-end data will have 0 for sorted_end_pairs and single_unmatched_pairs
        try:
            d["duplicate_rate"] = (
                d["duplicate_reads"] / ((d["sorted_end_pairs"] * 2) + (d["single_ends"] - d["single_unmatched_pairs"]))
            ) * 100.0
        except ZeroDivisionError:
            d["duplicate_rate"] = 0
            log.debug(f"Sambamba Markdup: zero division error for '{f['fn']}'")

        # Calculate some read counts from pairs - Paired End
        if d["sorted_end_pairs"] > 0:
            d["total_sorted_paired_end_reads"] = (d["sorted_end_pairs"] * 2) - d["duplicate_reads"]
            d["total_single_end_reads"] = d["single_ends"] - (2 * d["single_unmatched_pairs"])
            d["total_single_unmatched_reads"] = 2 * d["single_unmatched_pairs"]
        # Calculate some read counts from pairs - Single End
        else:
            d["total_sorted_paired_end_reads"] = 0
            d["total_single_end_reads"] = d["single_ends"] - d["duplicate_reads"]
            d["total_single_unmatched_reads"] = 0

        return d

    def markdup_general_stats_table(self):
        """
        Take parsed stats from sambamba markdup to general stats table at the top of the report.

        Append to the dictionaries from the top-level Sambamba module script, so that all columns
        are added in a single function call across all Sambamba submodules.
        """

        self.general_stats_headers["duplicate_rate"] = {
            "title": "% Dups",
            "description": "Rate of Duplication per Sample",
            "scale": "RdYlGn-rev",
            "format": "{:,.0f}",
            "suffix": "%",
        }
        for s_name, data in self.markdup_data.items():
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = {}
            self.general_stats_data[s_name]["duplicate_rate"] = data["duplicate_rate"]

    def markdup_section(self):
        """
        Add markdup statistics as bar graph to multiQC report.
        """

        # plot these categories, but not duplicate rate.
        cats = {
            "total_sorted_paired_end_reads": {"name": "Total Sorted Paired End Reads"},
            "total_single_end_reads": {"name": "Total Single End Reads"},
            "total_single_unmatched_reads": {"name": "Total Single Unmatched Reads"},
            "duplicate_reads": {"name": "Total Duplicate Reads"},
        }
        config = {
            "id": "SambambaMarkdupBargraph",
            "title": "Sambamba Markdup: Duplicate Counts",
            "ylab": "# Reads",
            "cpswitch_c_active": False,
        }

        self.add_section(
            name="Sambamba Markdup",
            anchor="SambambaMarkdup",
            description="Number of reads, categorised by duplication state. **Pair counts are doubled** - see help text for details.",
            helptext="""
                The table in the Picard metrics file contains some columns referring read pairs and some referring to single reads.

                To make the numbers in this plot sum correctly, values referring to pairs are doubled according to the scheme below:

                * `Total Sorted Paired End Reads = (Sorted End Pairs * 2) - Duplicate Reads`
                * `Total Single End Reads = Single Ends - (2 * Single Unmatched Pairs)`
                * `Total Single Unmatched Reads = 2 * Single Unmatched Pairs`
            """,
            plot=bargraph.plot(self.markdup_data, cats, config),
        )
