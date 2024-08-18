"""MultiQC module to parse output from ganon classify"""

import logging
import re
from typing import Union, Dict, Tuple, Optional

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module takes summary statistics from a file containing stdout from `ganon classify`.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Ganon",
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
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} reports")
        self.write_data_file(data_by_sample, "ganon")

        self.stats_table(data_by_sample)

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

    def stats_table(self, data_by_sample):
        """Ganon general stats table"""
        headers: Dict[str, Dict] = {
            "reads_processed": {
                "title": "Processed",
                "description": "Number of reads processed by ganon",
                "scale": "Greens",
                "shared_key": "read_count",
            },
            "mbp_processed": {
                "title": "Bases processed",
                "description": "Number of bases processed by ganon in Mbp",
                "scale": "Purples",
                "hidden": True,
            },
            "reads_classified": {
                "title": "Classified",
                "description": "Number of processed reads taxonomically classified",
                "scale": "Blues",
                "shared_key": "read_count",
                "hidden": True,
            },
            "reads_classified_pc": {
                "title": "Classified",
                "description": "Percent of processed reads taxonomically classified",
                "scale": "RdYlGn",
                "suffix": "%",
            },
            "unique_matches": {
                "title": "Unique matches",
                "description": "Number of reads with unique matches",
                "scale": "Blues",
                "shared_key": "read_count",
                "hidden": True,
            },
            "unique_matches_pc": {
                "title": "Unique matches",
                "description": "Percent of reads with unique matches",
                "scale": "RdYlGn",
                "suffix": "%",
            },
            "multiple_matches": {
                "title": "Multiple matches",
                "description": "Number of reads with multiple matches",
                "scale": "Blues",
                "shared_key": "read_count",
                "hidden": True,
            },
            "multiple_matches_pc": {
                "title": "Multiple matches",
                "description": "Percent of reads with multiple matches",
                "scale": "RdYlGn",
                "suffix": "%",
            },
            "overall_matches": {
                "title": "Overall matches",
                "description": "Number of matches overall",
                "scale": "Greens",
                "hidden": True,
            },
            "match_to_read": {
                "title": "Avg matches per read",
                "description": "Average number matches for each read",
                "scale": "Blues",
            },
            "reads_unclassified": {
                "title": "Unclassified",
                "description": "Number of reads processed reads without taxonomic classification",
                "scale": "Reds",
                "shared_key": "read_count",
                "hidden": True,
            },
            "reads_unclassified_pc": {
                "title": "Unclassified",
                "description": "Percent of processed reads without any taxonomic classification",
                "scale": "RdYlGn",
                "suffix": "%",
            },
            "taxonomic_entries_reported": {
                "title": "Taxa identified",
                "description": "Number of taxonomic entries in ganon report",
                "scale": "Greens",
            },
            "taxonomic_entries_removed_rank_filter": {
                "title": "Taxa rank removed",
                "description": "Number of taxonomic entries removed due to rank filter",
                "scale": "Purples",
                "hidden": True,
            },
            "taxonomic_entries_removed_mincount_filter": {
                "title": "Taxa min count removed",
                "description": "Number of taxonomic entries removed due to minimum count filter",
                "scale": "Reds",
                "hidden": True,
            },
        }

        self.add_section(
            name="Classify stats",
            anchor="ganon-stats-section",
            description="Summary of Ganon (classify) statistics",
            plot=table.plot(
                data_by_sample,
                headers,
                pconfig={
                    "id": "ganon-stats-table",
                    "title": "Ganon (classify) stats",
                },
            ),
        )

        general_stats_headers: Dict[str, Dict] = {}
        for k in [
            "reads_processed",
            "reads_classified_pc",
            "unique_matches_pc",
            "multiple_matches_pc",
            "match_to_read",
            "reads_unclassified_pc",
            "taxonomic_entries_reported",
        ]:
            general_stats_headers[k] = headers[k]
            general_stats_headers[k]["hidden"] = k not in [
                "taxonomic_entries_reported",
            ]

        self.general_stats_addcols(data_by_sample, general_stats_headers)

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
