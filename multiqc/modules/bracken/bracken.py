#!/usr/bin/env python

""" MultiQC module to parse output from bracken """

from __future__ import print_function
from collections import OrderedDict
import os
import logging
import re

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ Bracken module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Bracken",
            anchor="bracken",
            href="https://ccb.jhu.edu/software/bracken/",
            info="is a highly accurate statistical method that computes the abundance of species in DNA sequences from a metagenomics sample.",
        )

        self.t_ranks = OrderedDict()
        self.t_ranks["S"] = "Species"
        self.t_ranks["G"] = "Genus"
        self.t_ranks["F"] = "Family"
        self.t_ranks["O"] = "Order"
        self.t_ranks["C"] = "Class"
        self.t_ranks["P"] = "Phylum"
        self.t_ranks["K"] = "Kingdom"
        self.t_ranks["D"] = "Domain"
        self.t_ranks["R"] = "Root"
        # self.t_ranks['U'] = 'Unclassified'

        # Find and load any bracken reports
        self.bracken_raw_data = dict()
        for f in self.find_log_files("bracken", filehandles=True):
            self.parse_logs(f)

        # Filter to strip out ignored sample names
        self.bracken_raw_data = self.ignore_samples(self.bracken_raw_data)

        if len(self.bracken_raw_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.bracken_raw_data)))

        # Sum counts across all samples, so that we can pick top 5
        self.bracken_total_pct = dict()
        self.bracken_total_counts = dict()
        self.bracken_sample_total_readcounts = dict()
        self.sample_total_readcounts()
        self.sum_sample_counts()

        self.general_stats_cols()
        self.top_five_barplot()

    def parse_logs(self, f):
        """
        Parses a bracken report output file

        1. Percentage of fragments covered by the clade rooted at this taxon
        2. Number of fragments covered by the clade rooted at this taxon
        3. Number of fragments assigned directly to this taxon
        4. A rank code, indicating:
            * (U)nclassified
            * (R)oot
            * (D)omain
            * (K)ingdom
            * (P)hylum
            * (C)lass
            * (O)rder
            * (F)amily
            * (G)enus
            * (S)pecies
           Taxa that are not at any of these 10 ranks have a rank code that is
           formed by using the rank code of the closest ancestor rank with
           a number indicating the distance from that rank.  E.g., "G2" is a
           rank code indicating a taxon is between genus and species and the
           grandparent taxon is at the genus rank.
        5. NCBI taxonomic ID number
        6. Indented scientific name
        """

        # Search regexes for stats
        br_regex = re.compile(r"^(\d{1,3}\.\d{1,2})\t(\d+)\t(\d+)\t([\dRUDKPCOFGS-]{1,3})\t(\d+)(\s+)(.+)")
        data = []
        for l in f["f"]:
            match = br_regex.search(l)
            if match:
                row = {
                    "percent": float(match.group(1)),
                    "counts_rooted": int(match.group(2)),
                    "counts_direct": int(match.group(3)),
                    "rank_code": match.group(4),
                    "tax_id": int(match.group(5)),
                    "num_spaces": len(match.group(6)),
                    "classif": match.group(7),
                }
                data.append(row)

        self.bracken_raw_data[f["s_name"]] = data

    def sample_total_readcounts(self):
        """ Compute the total read counts for each sample """

        for s_name, data in self.bracken_raw_data.items():
            self.bracken_sample_total_readcounts[s_name] = 0
            for row in data:
                self.bracken_sample_total_readcounts[s_name] += \
                    row['counts_direct']

    def sum_sample_counts(self):
        """ Sum counts across all samples for bracken data """

        # Sum the percentages for each taxa across all samples
        # Allows us to pick top-5 for each rank
        # Use percentages instead of counts so that deeply-sequences samples
        # are not unfairly over-represented
        for s_name, data in self.bracken_raw_data.items():
            for row in data:

                # Convenience vars that are easier to read
                rank_code = row["rank_code"]
                classif = row["classif"]

                # Skip anything that doesn't exactly fit a tax rank level
                if row["rank_code"] == "-" or any(c.isdigit() for c in row["rank_code"]):
                    continue

                if rank_code not in self.bracken_total_pct:
                    self.bracken_total_pct[rank_code] = dict()
                    self.bracken_total_counts[rank_code] = dict()

                if classif not in self.bracken_total_pct[rank_code]:
                    self.bracken_total_pct[rank_code][classif] = 0
                    self.bracken_total_counts[rank_code][classif] = 0
                self.bracken_total_pct[rank_code][classif] += row["percent"]
                self.bracken_total_counts[rank_code][classif] += row["counts_rooted"]

    def general_stats_cols(self):
        """ Add a couple of columns to the General Statistics table """

        # Get top taxa in most specific taxa rank that we have
        top_five = []
        top_rank_code = None
        top_rank_name = None
        for rank_code, rank_name in self.t_ranks.items():
            try:
                sorted_pct = sorted(self.bracken_total_pct[rank_code].items(), key=lambda x: x[1], reverse=True)
                for classif, pct_sum in sorted_pct[:5]:
                    top_five.append(classif)
                top_rank_code = rank_code
                top_rank_name = rank_name
                break
            except KeyError:
                # No species-level data found etc
                pass

        top_one_hkey = "% {}".format(top_five[0])

        # Column headers
        headers = OrderedDict()
        headers[top_one_hkey] = {
            "title": top_one_hkey,
            "description": "Percentage of reads that were the top {} over all samples - {}".format(
                top_rank_name, top_five[0]
            ),
            "suffix": "%",
            "max": 100,
            "scale": "PuBuGn",
        }
        headers["% Top 5"] = {
            "title": "% Top 5 {}".format(top_rank_name),
            "description": "Percentage of reads that were classified by one of the top 5 {} ({})".format(
                top_rank_name, ", ".join(top_five)
            ),
            "suffix": "%",
            "max": 100,
            "scale": "PuBu",
        }
        headers["% Unclassified"] = {
            "title": "% Unclassified",
            "description": "Percentage of reads that were unclassified",
            "suffix": "%",
            "max": 100,
            "scale": "OrRd",
        }

        # Get table data
        tdata = {}
        for s_name, d in self.bracken_raw_data.items():
            tdata[s_name] = {}
            for row in d:
                percent = (row['counts_rooted'] / self.bracken_sample_total_readcounts[s_name]) * 100
                if row["rank_code"] == "U":
                    tdata[s_name]['% Unclassified'] = percent
                if row["rank_code"] == top_rank_code and row["classif"] in top_five:
                    tdata[s_name]['% Top 5'] = percent + tdata[s_name].get('% Top 5', 0)
                if row["rank_code"] == top_rank_code and row["classif"] == top_five[0]:
                    tdata[s_name][top_one_hkey] = percent

            if top_one_hkey not in tdata[s_name]:
                tdata[s_name][top_one_hkey] = 0

        self.general_stats_addcols(tdata, headers)

    def top_five_barplot(self):
        """ Add a bar plot showing the top-5 from each taxa rank """

        pd = []
        cats = list()
        pconfig = {
            "id": "bracken-topfive-plot",
            "title": "Bracken 2: Top taxa",
            "ylab": "Number of fragments",
            "data_labels": list(self.t_ranks.values()),
        }

        for rank_code, rank_name in self.t_ranks.items():
            rank_cats = OrderedDict()
            rank_data = dict()

            # Loop through the summed tax percentages to get the top 5 across all samples
            try:
                sorted_pct = sorted(self.bracken_total_pct[rank_code].items(), key=lambda x: x[1], reverse=True)
            except KeyError:
                # Taxa rank not found in this sample
                continue
            i = 0
            counts_shown = {}
            for classif, pct_sum in sorted_pct:
                i += 1
                if i > 5:
                    break
                rank_cats[classif] = {"name": classif}
                # Pull out counts for this rank + classif from each sample
                for s_name, d in self.bracken_raw_data.items():
                    if s_name not in rank_data:
                        rank_data[s_name] = dict()
                    if s_name not in counts_shown:
                        counts_shown[s_name] = 0
                    for row in d:
                        if row["rank_code"] == rank_code:
                            if row["classif"] == classif:
                                if classif not in rank_data[s_name]:
                                    rank_data[s_name][classif] = 0
                                rank_data[s_name][classif] += row["counts_rooted"]
                                counts_shown[s_name] += row["counts_rooted"]

            # Add in unclassified reads and "other" - we presume from other species etc.
            for s_name, d in self.bracken_raw_data.items():
                for row in d:
                    if row["rank_code"] == "U":
                        rank_data[s_name]["U"] = row["counts_rooted"]
                        counts_shown[s_name] += row["counts_rooted"]
                rank_data[s_name]["other"] = self.bracken_sample_total_readcounts[s_name] - counts_shown[s_name]

                # This should never happen... But it does sometimes if the total read count is a bit off
                if rank_data[s_name]["other"] < 0:
                    log.debug(
                        "Found negative 'other' count for {} ({}): {}".format(
                            s_name, self.t_ranks[rank_code], rank_data[s_name]["other"]
                        )
                    )
                    rank_data[s_name]["other"] = 0

            rank_cats["other"] = {"name": "Other", "color": "#cccccc"}
            rank_cats["U"] = {"name": "Unclassified", "color": "#d4949c"}

            cats.append(rank_cats)
            pd.append(rank_data)

        self.add_section(
            name="Top taxa",
            anchor="bracken-topfive",
            description="The number of reads falling into the top 5 taxa across different ranks.",
            helptext="""
                To make this plot, the percentage of each sample assigned to a given taxa is summed across all samples.
                The counts for these top five taxa are then plotted for each of the 9 different taxa ranks.
                The unclassified count is always shown across all taxa ranks.

                The total number of reads is approximated by dividing the number of `unclassified` reads by the percentage of
                the library that they account for.
                Note that this is only an approximation, and that bracken percentages don't always add to exactly 100%.

                The category _"Other"_ shows the difference between the above total read count and the sum of the read counts
                in the top 5 taxa shown + unclassified. This should cover all taxa _not_ in the top 5, +/- any rounding errors.

                Note that any taxon that does not exactly fit a taxon rank (eg. `-` or `G2`) is ignored.
            """,
            plot=bargraph.plot(pd, cats, pconfig),
        )
