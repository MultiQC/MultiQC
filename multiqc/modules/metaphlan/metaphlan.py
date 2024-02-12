""" MultiQC module to parse output from MetaPhlAn """


import logging
import re

from multiqc.utils import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """MetaPhlAn module"""

    def __init__(
        # Initialise the parent object
        self,
        name="MetaPhlAn",
        anchor="metaphlan",
        href="https://github.com/biobakery/MetaPhlAn",
        info="is a computational tool for profiling the composition of microbial communities from metagenomic shotgun sequencing data.",
        doi="10.1038/s41587-023-01688-w",
    ):
        super(MultiqcModule, self).__init__(
            name=name,
            anchor=anchor,
            href=href,
            info=info,
            doi=doi,
        )
        # Custom options from user config that can overwrite base module values
        self.t_ranks = {
            "s": "Species",
            "g": "Genus",
            "f": "Family",
            "o": "Order",
            "c": "Class",
            "p": "Phylum",
            "k": "Kingdom",
        }

        self.top_n = getattr(config, "metaphlan", {}).get("top_n", 10)

        self.metaphlan_raw_data = dict()
        for f in self.find_log_files("metaphlan", filehandles=True):
            f["f"].seek(0)
            self.parse_logs(f)
            self.add_data_source(f)

        self.metaphlan_raw_data = self.ignore_samples(self.metaphlan_raw_data)

        if len(self.metaphlan_raw_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.metaphlan_raw_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.write_data_file(self.metaphlan_raw_data, f"multiqc_{self.anchor}")

        # Sum percentages across all samples, so that we can pick top species
        self.metaphlan_total_pct = dict()
        self.sum_sample_counts()

        self.general_stats_cols()
        self.top_taxa_barplot()

    def parse_logs(self, f):
        """
        Parses a metaphlan output .txt file to extract the following information

        1. A rank code, indicating:
            * (k)ingdom
            * (p)hylum
            * (c)lass
            * (o)rder
            * (f)amily
            * (g)enus
            * (s)pecies
        2. Clade name
        3. Relative abundance of this clade
        """

        # Search regexes for stats
        regex_kingdom_level = re.compile(r"([kpcofgs])\_\_(\w+)\t(\d+)\t(\d{1,3}\.\d{1,5})")
        regex_rel_abundance = re.compile(r"\t(\d{1,3}\.\d{1,5})")
        regex_last_taxid = re.compile(r"\|(\d+)\t")
        regex_last_tax_rank = re.compile(r"\|([kpcofgs])\_\_(\w+)\t")
        regex_tax_unclass = re.compile(r"\|([kpcofgs])\_\_(\w+)\_(unclassified)\t")
        regex_comments = re.compile(r"\#")
        data = []
        for line in f["f"]:
            match = regex_last_tax_rank.search(line)
            match_last_id = regex_last_taxid.search(line)
            match_per = regex_rel_abundance.search(line)
            match_first = regex_kingdom_level.search(line)
            match_unclassified = regex_tax_unclass.search(line)
            match_comments = regex_comments.search(line)
            if match_first:
                # Matches the first row (kingdom)
                row = {
                    "tax_rank": match_first.group(1),
                    "taxonomy": match_first.group(2),
                    "rel_abundance": float(match_first.group(4)),
                }
                data.append(row)
            elif match_unclassified and match_per:
                # Include unclassified reads on each rank
                row = {
                    "tax_rank": match_unclassified.group(1),
                    "taxonomy": match_unclassified.group(3),
                    "rel_abundance": float(match_per.group(1)),
                }
                data.append(row)
            elif match and match_last_id and match_per:
                # Match all other clades
                row = {
                    "tax_rank": match.group(1),
                    "taxonomy": match.group(2),
                    "rel_abundance": float(match_per.group(1)),
                }
                data.append(row)
            elif match_comments:
                continue
            else:
                log.debug(f"{f['s_name']}: Could not parse line: {line}")
        self.metaphlan_raw_data[f["s_name"]] = data

    def sum_sample_counts(self):
        """Sum relative abundance across all samples for metaphlan data"""

        # Sum the percentages for each taxa across all samples
        # Allows us to pick the top taxa for each rank
        for s_name, data in self.metaphlan_raw_data.items():
            for row in data:
                tax_rank = row["tax_rank"]
                taxonomy = row["taxonomy"]

                if tax_rank not in self.metaphlan_total_pct:
                    self.metaphlan_total_pct[tax_rank] = dict()

                if taxonomy not in self.metaphlan_total_pct[tax_rank]:
                    self.metaphlan_total_pct[tax_rank][taxonomy] = 0
                self.metaphlan_total_pct[tax_rank][taxonomy] += row["rel_abundance"]

    def general_stats_cols(self):
        """Add a couple of columns to the General Statistics table"""

        # Get top taxa in most specific taxa rank that we have
        top_taxa = []
        top_rank_code = None
        top_rank_name = None
        for rank_code, rank_name in self.t_ranks.items():
            try:
                sorted_pct = sorted(
                    self.metaphlan_total_pct[rank_code].items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
                for taxonomy, pct_sum in sorted_pct[: self.top_n]:
                    top_taxa.append(taxonomy)
                top_rank_code = rank_code
                top_rank_name = rank_name
                break
            except KeyError:
                # No species-level data found etc
                pass

        # Column headers
        headers = dict()

        top_one_hkey = None

        top_one_hkey = "% {}".format(top_taxa[0])
        headers[top_one_hkey] = {
            "title": top_one_hkey,
            "description": "Percentage of reads that were the top {} over all samples - {}".format(
                top_rank_name, top_taxa[0]
            ),
            "suffix": "%",
            "max": 100,
            "scale": "PuBuGn",
        }
        headers["% Top"] = {
            "title": f"% Top {self.top_n} {top_rank_name}",
            "description": f"Percentage of reads that were classified by one of the top-{self.top_n} {top_rank_name} ({', '.join(top_taxa)})",
            "suffix": "%",
            "max": 100,
            "scale": "PuBu",
        }
        # Get table data
        tdata = {}
        for s_name, d in self.metaphlan_raw_data.items():
            tdata[s_name] = {}
            for row in d:
                percent = row["rel_abundance"]
                if row["tax_rank"] == top_rank_code and row["taxonomy"] in top_taxa:
                    tdata[s_name]["% Top"] = percent + tdata[s_name].get("% Top", 0)
                if row["tax_rank"] == top_rank_code and row["taxonomy"] == top_taxa[0]:
                    tdata[s_name][top_one_hkey] = percent

            if top_one_hkey is not None and top_one_hkey not in tdata[s_name]:
                tdata[s_name][top_one_hkey] = 0

        self.general_stats_addcols(tdata, headers)

    def top_taxa_barplot(self):
        """Add a bar plot showing the top-N from each taxa rank"""

        pd = []
        cats = list()
        # Keeping track of encountered codes to display only tabs with available data
        found_rank_codes = set()

        for rank_code in self.t_ranks:
            rank_cats = dict()
            rank_data = dict()

            # Loop through the summed tax percentages to get the top-N across all samples
            try:
                sorted_pct = sorted(
                    self.metaphlan_total_pct[rank_code].items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
            except KeyError:
                # Taxa rank not found in this sample
                continue
            i = 0
            counts_shown = {}
            for taxonomy, pct_sum in sorted_pct:
                # Adding unclassified for each rank:
                if taxonomy == "unclassified":
                    rank_cats[taxonomy] = {"name": taxonomy}
                    for s_name, d in self.metaphlan_raw_data.items():
                        if s_name not in rank_data:
                            rank_data[s_name] = dict()
                        if s_name not in counts_shown:
                            counts_shown[s_name] = 0
                        for row in d:
                            if row["tax_rank"] == rank_code:
                                found_rank_codes.add(rank_code)
                                if row["taxonomy"] == taxonomy:
                                    if taxonomy not in rank_data[s_name]:
                                        rank_data[s_name][taxonomy] = 0
                                    rank_data[s_name][taxonomy] += row["rel_abundance"]
                                    counts_shown[s_name] += row["rel_abundance"]
                                    rank_data[s_name]["other"] = 100 - counts_shown[s_name]
                else:
                    # Add top n taxa for each rank
                    i += 1
                    if i > self.top_n:
                        # After top 5, keep looping to sum up unclassified
                        continue
                    rank_cats[taxonomy] = {"name": taxonomy}
                    # Pull out counts for this rank + classif from each sample
                    for s_name, d in self.metaphlan_raw_data.items():
                        if s_name not in rank_data:
                            rank_data[s_name] = dict()
                        if s_name not in counts_shown:
                            counts_shown[s_name] = 0

                        for row in d:
                            if row["tax_rank"] == rank_code:
                                found_rank_codes.add(rank_code)
                                if row["taxonomy"] == taxonomy:
                                    if taxonomy not in rank_data[s_name]:
                                        rank_data[s_name][taxonomy] = 0
                                    rank_data[s_name][taxonomy] += row["rel_abundance"]
                                    counts_shown[s_name] += row["rel_abundance"]
                                    rank_data[s_name]["other"] = 100 - counts_shown[s_name]
            # Add in unclassified reads and "other" - we presume from other species etc.
            for s_name, d in self.metaphlan_raw_data.items():
                # In case none of the top_n were in some sample:
                if s_name not in rank_data:
                    rank_data[s_name] = dict()
                if s_name not in counts_shown:
                    counts_shown[s_name] = 0
                rank_data[s_name]["other"] = 100 - counts_shown[s_name]

                # This should never happen... But it does sometimes if the total read count is a bit off
                if rank_data[s_name]["other"] < 0:
                    log.debug(
                        "Found negative 'other' count for {} ({}): {}".format(
                            s_name, self.t_ranks[rank_code], rank_data[s_name]["other"]
                        )
                    )
                    rank_data[s_name]["other"] = 0

            rank_cats["other"] = {"name": "Other", "color": "#cccccc"}
            rank_cats["unclassified"] = {"name": "Unclassified", "color": "#d4949c"}

            cats.append(rank_cats)
            pd.append(rank_data)

        pconfig = {
            "id": f"{self.anchor}-top-n-plot",
            "title": f"{self.name}: Top taxa",
            "ylab": "Number of fragments",
            "data_labels": [v for k, v in self.t_ranks.items() if k in found_rank_codes],
            "cpswitch": False,
        }

        self.add_section(
            name="Top taxa",
            anchor=f"{self.anchor}-top-n",
            description=f"The relative abundance of reads falling into the top {self.top_n} taxa across different ranks.",
            helptext=f"""
                To make this plot, the percentage of each sample assigned to a given taxa is summed across all samples.
                The counts for these top {self.top_n} taxa are then plotted for each of the 7 different taxa ranks.
                The unclassified percentages are shown per taxa rank.

                The category _"Other"_ shows the difference between 100% and the sum of the percent
                in the top {self.top_n} taxa shown + the percent that was unclassified. This should cover all taxa _not_ in the top {self.top_n}, +/- any rounding errors.

            """,
            plot=bargraph.plot(pd, cats, pconfig),
        )
