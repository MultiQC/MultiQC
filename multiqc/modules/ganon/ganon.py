""" MultiQC module to parse output from ganon classify """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, beeswarm
from multiqc.utils import config

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="ganon",
            anchor="ganon",
            href="https://pirovc.github.io/ganon/",
            info="ganon is developed for, but not limited, to the metagenomics classification problem: quickly assign sequence fragments to their closest reference among thousands of references via  Interleaved Bloom Filters of k-mer/minimizers.",
            doi="10.1093/bioinformatics/btaa458",
        )

        self.ganon_data = dict()

        for f in self.find_log_files("ganon"):
            self.parse_logs(f)

        self.calculate_entry_remainder()

        self.ganon_data = self.ignore_samples(self.ganon_data)

        if len(self.ganon_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.ganon_data)))
        print(self.ganon_data)

        self.write_data_file(self.ganon_data, "ganon")
        self.ganon_general_stats()

        self.barplot_reads_classified()
        self.barplot_reads_match_type()
        self.beeswarm_average_matches()
        self.barplot_taxonomic_entries()

    def parse_logs(self, f):
        for l in f["f"].splitlines():
            if l.startswith("--output-prefix"):
                ## find and set name - we don't clean as can never take from file name
                s_name = l.split()[1]

                ## check for duplicates
                if s_name in self.ganon_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.ganon_data[s_name] = {}

                self.add_data_source(f, s_name=s_name)

            if l.startswith("ganon-classify processed"):
                self.ganon_data[s_name]["reads_processed"] = int(l.split()[2])
                self.ganon_data[s_name]["mbp_processed"] = l.split()[4].lstrip("(").rstrip(")")

            if "reads classified" in l:
                self.ganon_data[s_name]["reads_classified"] = int(l.split()[1])
                self.ganon_data[s_name]["reads_classified_pc"] = l.split()[4].lstrip("(").rstrip(")").rstrip("%")

            if "with unique matches" in l:
                self.ganon_data[s_name]["unique_matches"] = int(l.split()[1])
                self.ganon_data[s_name]["unique_matches_pc"] = l.split()[5].lstrip("(").rstrip(")").rstrip("%")

            if "with multiple matches" in l:
                self.ganon_data[s_name]["multiple_matches"] = int(l.split()[1])
                self.ganon_data[s_name]["multiple_matches_pc"] = l.split()[5].lstrip("(").rstrip(")").rstrip("%")

            if "matches (avg" in l:
                self.ganon_data[s_name]["overall_matches"] = int(l.split()[1])
                self.ganon_data[s_name]["match_to_read"] = l.split()[4].lstrip("(").rstrip(")")

            if "reads unclassified" in l:
                self.ganon_data[s_name]["reads_unclassified"] = int(l.split()[1])
                self.ganon_data[s_name]["reads_unclassified_pc"] = l.split()[4].lstrip("(").rstrip(")").rstrip("%")

            if "entries reported" in l:
                self.ganon_data[s_name]["taxonomic_entries_reported"] = int(l.split()[1])

            if "removed not in --ranks" in l:
                self.ganon_data[s_name]["taxonomic_entries_removed_rank_filter"] = int(l.split()[1])

            if "removed with --min-count" in l:
                self.ganon_data[s_name]["taxonomic_entries_removed_mincount_filter"] = int(l.split()[1])

        return

    def calculate_entry_remainder(self):
        for s_name in self.ganon_data:
            self.ganon_data[s_name]["taxonomic_entries_retained"] = self.ganon_data[s_name][
                "taxonomic_entries_reported"
            ] - (
                self.ganon_data[s_name]["taxonomic_entries_removed_rank_filter"]
                + self.ganon_data[s_name]["taxonomic_entries_removed_mincount_filter"]
            )

    def ganon_general_stats(self):
        """ganon General Stats Table"""
        headers = OrderedDict()

        headers["reads_processed"] = {
            "title": "Reads Processed ({})".format(config.read_count_prefix),
            "description": "Number of reads processed by ganon ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }

        headers["mbp_processed"] = {
            "title": "Reads Processed (Mbp)",
            "description": "Number of reads processed by ganon in Mbp",
            "scale": "Purples",
        }

        headers["reads_classified"] = {
            "title": "Reads Classified ({})".format(config.read_count_prefix),
            "description": "Number of processed reads taxonomically classified ({})".format(config.read_count_prefix),
            "scale": "Blues",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }

        headers["reads_classified_pc"] = {
            "title": "Reads Classified Percent",
            "description": "Percent of processed reads taxonomically classified",
            "scale": "RdYlGn",
            "suffix": "%",
            "min": 0,
            "max": 100,
        }

        headers["unique_matches"] = {
            "title": "Unique Matches ({})".format(config.read_count_prefix),
            "description": "Number of reads with unique matches ({})".format(config.read_count_prefix),
            "scale": "Blues",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }

        headers["unique_matches_pc"] = {
            "title": "Unique Matches Percent",
            "description": "Percent of reads with unique matches",
            "scale": "RdYlGn",
            "suffix": "%",
            "min": 0,
            "max": 100,
        }

        headers["multiple_matches"] = {
            "title": "Multiple Matches ({})".format(config.read_count_prefix),
            "description": "Number of reads with multiple matches ({})".format(config.read_count_prefix),
            "scale": "Blues",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }

        headers["multiple_matches_pc"] = {
            "title": "Multiple Matches Percent",
            "description": "Percent of reads with multiple matches",
            "scale": "RdYlGn",
            "suffix": "%",
            "min": 0,
            "max": 100,
        }

        headers["overall_matches"] = {
            "title": "Overall Matches",
            "description": "Number of matches overall",
            "scale": "Greens",
        }

        headers["match_to_read"] = {
            "title": "Avg Matches Per Read",
            "description": "Average number matches for each read",
            "scale": "Blues",
        }

        headers["reads_unclassified"] = {
            "title": "Reads Unclassified ({})".format(config.read_count_prefix),
            "description": "Number of reads processed reads without taxonomic classification ({})".format(
                config.read_count_prefix
            ),
            "scale": "Blues",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }

        headers["reads_unclassified_pc"] = {
            "title": "Reads Unclassified Percent",
            "description": "Percent of processed reads without any taxonomic classification",
            "scale": "RdYlGn",
            "suffix": "%",
            "min": 0,
            "max": 100,
        }

        ## I think it is relatively unlikely to get millions of taxonomic entries,
        ## so not using multiplier here
        headers["taxonomic_entries_reported"] = {
            "title": "Nr. Taxa Identified",
            "description": "Number of taxonomic entries in ganon report",
            "scale": "Greens",
        }

        headers["taxonomic_entries_removed_rank_filter"] = {
            "title": "Nr. Taxa Rank Removed",
            "description": "Number of taxonomic entries removed due to rank filter",
            "scale": "Purples",
        }

        headers["taxonomic_entries_removed_mincount_filter"] = {
            "title": "Nr. Taxa Min Count Removed",
            "description": "Number of taxonomic entries removed due to minimum count filter",
            "scale": "Reds",
        }

        self.general_stats_addcols(self.ganon_data, headers)

    def barplot_reads_classified(self):
        """Barplot of total number of reads classified"""
        cats = OrderedDict()
        cats["reads_classified"] = {"name": "Classified reads", "color": "#7cb5ec"}
        cats["reads_unclassified"] = {"name": "Unclassified reads", "color": "#f7a35c"}
        config = {
            "id": "ganon-reads-classified",
            "title": "Ganon (classify): classified reads summary",
            "ylab": "Read Counts",
        }
        self.add_section(
            name="Reads classified",
            anchor="ganon-reads-classified",
            description="Summary of whether reads were able to be taxonomically assigned. Total should match number of reads processed.",
            plot=bargraph.plot(self.ganon_data, cats, config),
        )

    def barplot_reads_match_type(self):
        """Barplot of total number of reads classified"""
        cats = OrderedDict()
        cats["unique_matches"] = {"name": "Reads with unique matches", "color": "#7cb5ec"}
        cats["multiple_matches"] = {"name": "Reads with multiple matches", "color": "#f7a35c"}
        config = {
            "id": "ganon-reads-match-type",
            "title": "Ganon (classify): match type summary",
            "ylab": "Read Counts",
        }
        self.add_section(
            name="Match type distribution",
            anchor="ganon-reads-match-type",
            description="Summary of whether read hits were unique or had multiple matches. Total should match number of reads classified.",
            plot=bargraph.plot(self.ganon_data, cats, config),
        )

    def beeswarm_average_matches(self):
        """Barplot of total number of reads classified"""
        cats = OrderedDict()
        cats["match_to_read"] = {"name": "Average match to read", "color": "#7cb5ec"}
        config = {
            "id": "ganon-reads-match-to-read-ratio",
            "title": "Ganon (classify): average number of matches per read",
            "ylab": "Average reads per match",
        }
        self.add_section(
            name="Average match to read",
            anchor="ganon-reads-match-type",
            description="Summary of how many taxonomic matches each read had on average.",
            plot=beeswarm.plot(self.ganon_data, cats, config),
        )

    def barplot_taxonomic_entries(self):
        """Barplot of total number of reads classified"""
        cats = OrderedDict()
        cats["taxonomic_entries_retained"] = {"name": "Retained Taxonomic assignments", "color": "#7cb5ec"}
        cats["taxonomic_entries_removed_rank_filter"] = {"name": "Rank filter removed", "color": "#f7a35c"}
        cats["taxonomic_entries_removed_mincount_filter"] = {"name": "Min. count filter removed", "color": "#fb9a99"}

        config = {
            "id": "ganon-taxonomic-entries",
            "title": "Ganon (classify): distribution of taxonomic entries through filtering",
            "ylab": "Entries",
        }
        self.add_section(
            name="Taxonomic entries",
            anchor="ganon-taxonomic-entries",
            description="Summary of how many taxa were identified overall and removed through filtering. Total should match Nr. Taxa Identified.",
            plot=bargraph.plot(self.ganon_data, cats, config),
        )
