"""MultiQC module to parse output from ganon classify"""

import logging
import re
from collections import OrderedDict
from typing import Union, Dict, Tuple, Optional

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, violin
from multiqc import config

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module takes summary statistics from a file containing stdout from `ganon classify`.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ganon",
            anchor="ganon",
            href="https://pirovc.github.io/ganon/",
            info="Metagenomics classification: quickly assigns sequence fragments to their closest reference among thousands of references via Interleaved Bloom Filters of k-mer/minimizers.",
            doi="10.1093/bioinformatics/btaa458",
        )

        data_by_sample: Dict[str, Dict[str, Union[str, int, float, None]]] = dict()
        for f in self.find_log_files("ganon"):
            s_name, data = self.parse_log(f)
            if not s_name:
                continue
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            data_by_sample[s_name] = data

        self.calculate_entry_remainder(data_by_sample)

        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise UserWarning
        log.info(f"Found {len(data_by_sample)} reports")
        self.write_data_file(data_by_sample, "ganon")

        self.ganon_general_stats(data_by_sample)

        self.barplot_reads_classified(data_by_sample)
        self.barplot_reads_match_type(data_by_sample)
        self.violin_average_matches(data_by_sample)
        self.barplot_taxonomic_entries(data_by_sample)

    def parse_log(self, f) -> Tuple[Optional[str], Dict[str, Union[str, int, float, None]]]:
        s_name = None
        version = None
        data: Dict[str, Union[str, int, float, None]] = dict()
        for line in f["f"].splitlines():
            if s_name is None:  # We are somewhere in the file header
                if version is None:
                    m = re.search(r"^\s+_\|\s+v. (\d\.\d.\d)", line)
                    if m:
                        version = m
                        continue

                if line.startswith("--output-prefix"):
                    if s_name is not None:
                        log.debug(f"Duplicate sample name found within the same file {f['fn']}! Overwriting: {s_name}")
                    s_name = self.clean_s_name(line.split()[1])
                    self.add_data_source(f, s_name=s_name)
                    continue

            if line.startswith("ganon-classify processed"):
                data["reads_processed"] = int(line.split()[2])
                data["mbp_processed"] = line.split()[4].lstrip("(").rstrip(")")

            elif "reads classified" in line:
                data["reads_classified"] = int(line.split()[1])
                data["reads_classified_pc"] = line.split()[4].lstrip("(").rstrip(")").rstrip("%")

            elif "with unique matches" in line:
                data["unique_matches"] = int(line.split()[1])
                data["unique_matches_pc"] = line.split()[5].lstrip("(").rstrip(")").rstrip("%")

            elif "with multiple matches" in line:
                data["multiple_matches"] = int(line.split()[1])
                data["multiple_matches_pc"] = line.split()[5].lstrip("(").rstrip(")").rstrip("%")

            elif "matches (avg" in line:
                data["overall_matches"] = int(line.split()[1])
                data["match_to_read"] = line.split()[4].lstrip("(").rstrip(")")

            elif "reads unclassified" in line:
                data["reads_unclassified"] = int(line.split()[1])
                data["reads_unclassified_pc"] = line.split()[4].lstrip("(").rstrip(")").rstrip("%")

            elif "entries reported" in line:
                data["taxonomic_entries_reported"] = int(line.split()[1])

            elif "removed not in --ranks" in line:
                data["taxonomic_entries_removed_rank_filter"] = int(line.split()[1])

            elif "removed with --min-count" in line:
                data["taxonomic_entries_removed_mincount_filter"] = int(line.split()[1])

        if version is not None:
            self.add_software_version(version.group(1), sample=s_name)

        return s_name, data

    @staticmethod
    def calculate_entry_remainder(data_by_sample):
        for s_name in data_by_sample:
            data_by_sample[s_name]["taxonomic_entries_retained"] = data_by_sample[s_name][
                "taxonomic_entries_reported"
            ] - (
                data_by_sample[s_name]["taxonomic_entries_removed_rank_filter"]
                + data_by_sample[s_name]["taxonomic_entries_removed_mincount_filter"]
            )

    def ganon_general_stats(self, data_by_sample):
        """Ganon general stats table"""
        headers = {
            "reads_processed": {
                "title": "Reads Processed ({})".format(config.read_count_prefix),
                "description": "Number of reads processed by ganon ({})".format(config.read_count_prefix),
                "scale": "Greens",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "mbp_processed": {
                "title": "Reads Processed (Mbp)",
                "description": "Number of reads processed by ganon in Mbp",
                "scale": "Purples",
            },
            "reads_classified": {
                "title": "Reads Classified ({})".format(config.read_count_prefix),
                "description": "Number of processed reads taxonomically classified ({})".format(
                    config.read_count_prefix
                ),
                "scale": "Blues",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "reads_classified_pc": {
                "title": "Reads Classified Percent",
                "description": "Percent of processed reads taxonomically classified",
                "scale": "RdYlGn",
                "suffix": "%",
                "min": 0,
                "max": 100,
            },
            "unique_matches": {
                "title": "Unique Matches ({})".format(config.read_count_prefix),
                "description": "Number of reads with unique matches ({})".format(config.read_count_prefix),
                "scale": "Blues",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "unique_matches_pc": {
                "title": "Unique Matches Percent",
                "description": "Percent of reads with unique matches",
                "scale": "RdYlGn",
                "suffix": "%",
                "min": 0,
                "max": 100,
            },
            "multiple_matches": {
                "title": "Multiple Matches ({})".format(config.read_count_prefix),
                "description": "Number of reads with multiple matches ({})".format(config.read_count_prefix),
                "scale": "Blues",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "multiple_matches_pc": {
                "title": "Multiple Matches Percent",
                "description": "Percent of reads with multiple matches",
                "scale": "RdYlGn",
                "suffix": "%",
                "min": 0,
                "max": 100,
            },
            "overall_matches": {
                "title": "Overall Matches",
                "description": "Number of matches overall",
                "scale": "Greens",
            },
            "match_to_read": {
                "title": "Avg Matches Per Read",
                "description": "Average number matches for each read",
                "scale": "Blues",
            },
            "reads_unclassified": {
                "title": "Reads Unclassified ({})".format(config.read_count_prefix),
                "description": "Number of reads processed reads without taxonomic classification ({})".format(
                    config.read_count_prefix
                ),
                "scale": "Reds",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "reads_unclassified_pc": {
                "title": "Reads Unclassified Percent",
                "description": "Percent of processed reads without any taxonomic classification",
                "scale": "RdYlGn",
                "suffix": "%",
                "min": 0,
                "max": 100,
            },
            # I think it is relatively unlikely to get millions of taxonomic entries,
            # so not using multiplier here
            "taxonomic_entries_reported": {
                "title": "Nr. Taxa Identified",
                "description": "Number of taxonomic entries in ganon report",
                "scale": "Greens",
            },
            "taxonomic_entries_removed_rank_filter": {
                "title": "Nr. Taxa Rank Removed",
                "description": "Number of taxonomic entries removed due to rank filter",
                "scale": "Purples",
            },
            "taxonomic_entries_removed_mincount_filter": {
                "title": "Nr. Taxa Min Count Removed",
                "description": "Number of taxonomic entries removed due to minimum count filter",
                "scale": "Reds",
            },
        }

        self.general_stats_addcols(data_by_sample, headers)

    def barplot_reads_classified(self, data_by_sample):
        """Barplot of total number of reads classified"""
        self.add_section(
            name="Reads classified",
            anchor="ganon-reads-classified-section",
            description="Summary of whether reads were able to be taxonomically assigned. Total should match number of reads processed.",
            plot=bargraph.plot(
                data_by_sample,
                cats={
                    "reads_classified": {"name": "Classified reads", "color": "#7cb5ec"},
                    "reads_unclassified": {"name": "Unclassified reads", "color": "#f7a35c"},
                },
                pconfig={
                    "id": "ganon-reads-classified-plot",
                    "title": "Ganon (classify): classified reads summary",
                    "ylab": "Read Counts",
                },
            ),
        )

    def barplot_reads_match_type(self, data_by_sample):
        """Barplot of total number of reads classified"""
        self.add_section(
            name="Match type distribution",
            anchor="ganon-reads-match-type-section",
            description="Summary of whether read hits were unique or had multiple matches. Total should match number of reads classified.",
            plot=bargraph.plot(
                data_by_sample,
                cats={
                    "unique_matches": {
                        "name": "Reads with unique matches",
                        "color": "#7cb5ec",
                    },
                    "multiple_matches": {
                        "name": "Reads with multiple matches",
                        "color": "#f7a35c",
                    },
                },
                pconfig={
                    "id": "ganon-reads-match-type-plot",
                    "title": "Ganon (classify): match type summary",
                    "ylab": "Read Counts",
                },
            ),
        )

    def violin_average_matches(self, data_by_sample):
        """Barplot of total number of reads classified"""
        self.add_section(
            name="Average match to read",
            anchor="ganon-reads-match-to-read-ratio-section",
            description="Summary of how many taxonomic matches each read had on average.",
            plot=bargraph.plot(
                data_by_sample,
                cats={"match_to_read": {"name": "Average match to read", "color": "#7cb5ec"}},
                pconfig={
                    "id": "ganon-reads-match-to-read-ratio-plot",
                    "title": "Ganon (classify): average number of matches per read",
                    "ylab": "Average reads per match",
                },
            ),
        )

    def barplot_taxonomic_entries(self, data_by_sample):
        """Barplot of total number of reads classified"""
        self.add_section(
            name="Taxonomic entries",
            anchor="ganon-taxonomic-entries-section",
            description="Summary of how many taxa were identified overall and removed through filtering. Total should match Nr. Taxa Identified.",
            plot=bargraph.plot(
                data_by_sample,
                cats={
                    "taxonomic_entries_retained": {
                        "name": "Retained Taxonomic assignments",
                        "color": "#7cb5ec",
                    },
                    "taxonomic_entries_removed_rank_filter": {
                        "name": "Rank filter removed",
                        "color": "#f7a35c",
                    },
                    "taxonomic_entries_removed_mincount_filter": {
                        "name": "Min. count filter removed",
                        "color": "#fb9a99",
                    },
                },
                pconfig={
                    "id": "ganon-taxonomic-entries-plot",
                    "title": "Ganon (classify): distribution of taxonomic entries through filtering",
                    "ylab": "Entries",
                },
            ),
        )
