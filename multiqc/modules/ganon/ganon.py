""" MultiQC module to parse output from ganon classify """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph
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

        ## Define the main ganon multiqc data object
        self.ganon_data = dict()

        for f in self.find_log_files("ganon"):
            self.parse_logs(f)

        # self.ganon_data = self.ignore_samples(self.ganon_data)

        # if len(self.ganon_data) == 0:
        #     raise UserWarning

        log.info("Found {} reports".format(len(self.ganon_data)))
        print(self.ganon_data)

        self.write_data_file(self.ganon_data, "ganon")
        self.ganon_general_stats()

        ## TODO: Write to file and printing functions

        # self.ganon_beeswarm_plot()

    def parse_logs(self, f):
        for l in f["f"].splitlines():
            if l.startswith("--output-prefix"):
                ## find and set name
                s_name = l.split()[1]

                ## check for duplicates
                if s_name in self.ganon_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.ganon_data[s_name] = {}

                ## log the file name against the sample name(?)
                ##self.add_data_source(f, s_name=s_name)

            ## TODO Functionise?

            if l.startswith("ganon-classify processed"):
                self.ganon_data[s_name]["reads_processed"] = l.split()[2]
                self.ganon_data[s_name]["mbp_processed"] = l.split()[4].lstrip("(").rstrip(")")

            if "reads classified" in l:
                self.ganon_data[s_name]["reads_classified"] = l.split()[1]
                self.ganon_data[s_name]["reads_classified_pc"] = l.split()[4].lstrip("(").rstrip(")").rstrip("%")

            if "with unique matches" in l:
                self.ganon_data[s_name]["unique_matches"] = l.split()[1]
                self.ganon_data[s_name]["unique_matches_pc"] = l.split()[5].lstrip("(").rstrip(")").rstrip("%")

            if "with multiple matches" in l:
                self.ganon_data[s_name]["multiple_matches"] = l.split()[1]
                self.ganon_data[s_name]["multiple_matches_pc"] = l.split()[5].lstrip("(").rstrip(")").rstrip("%")

            if "matches (avg" in l:
                self.ganon_data[s_name]["overall_matches"] = l.split()[1]
                self.ganon_data[s_name]["match_to_read"] = l.split()[4].lstrip("(").rstrip(")")

            if "reads unclassified" in l:
                self.ganon_data[s_name]["reads_unclassified"] = l.split()[1]
                self.ganon_data[s_name]["reads_unclassified_pc"] = l.split()[4].lstrip("(").rstrip(")").rstrip("%")

            if "entries reported" in l:
                self.ganon_data[s_name]["taxonomic_entries_reported"] = l.split()[1]

            if "removed not in --ranks" in l:
                self.ganon_data[s_name]["taxonomic_entries_removed_rank_filter"] = l.split()[1]

            if "removed with --min-count" in l:
                self.ganon_data[s_name]["taxonomic_entries_removed_mincount_filter"] = l.split()[1]

        return

    def ganon_general_stats(self):
        """ganon General Stats Table"""
        headers = OrderedDict()

        ## TODO Fix read multiplier see error message, something to do with float
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


## TODO: barplot of reads classified + unclassified (total matches reads processed)
## TODO: barplot of reads with unique and multiple matches (total matches reads classified)
## TODO: barplot: total matches
## TODO: barplot: average matches to reads
## TODO: barplot: filter removed + filter removed minus total entries (remainder == confident etnries)

## TODO: to make testdataset submission
