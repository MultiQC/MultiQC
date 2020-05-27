#!/usr/bin/env python

""" MultiQC module to parse output from kraken """

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
    """ Kraken module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name = 'Kraken',
            anchor = 'kraken',
            href = "https://ccb.jhu.edu/software/kraken/",
            info = "is a taxonomic classification tool that uses exact k-mer matches to find the lowest common ancestor (LCA) of a given sequence."
        )

        # Find and load any kraken reports
        self.kraken_raw_data = dict()
        for f in self.find_log_files('kraken', filehandles=True):
            self.parse_logs(f)

        # Filter to strip out ignored sample names
        self.kraken_raw_data = self.ignore_samples(self.kraken_raw_data)

        if len(self.kraken_raw_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.kraken_raw_data)))

        # Sum counts across all samples, so that we can pick top 5
        self.kraken_total_pct = dict()
        self.kraken_total_counts = dict()
        self.kraken_sample_total_readcounts = dict()
        self.sum_sample_counts()

        self.top_five_barplot()


    def parse_logs(self, f):
        """
        Parses a kraken report output file

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
        k2_regex = re.compile(r"^\s{1,2}(\d{1,2}\.\d{1,2})\t(\d+)\t(\d+)\t([\dUDKPCOFGS-]{1,3})\t(\d+)(\s+)(.+)")
        data = []
        for l in f['f']:
            match = k2_regex.search(l)
            if match:
                row = {
                    'percent': float(match.group(1)),
                    'counts_rooted': int(match.group(2)),
                    'counts_direct': int(match.group(3)),
                    'rank_code': match.group(4),
                    'tax_id': int(match.group(5)),
                    'num_spaces': len(match.group(6)),
                    'classif': match.group(7)
                }
                data.append(row)

        self.kraken_raw_data[f['s_name']] = data

    def sum_sample_counts(self):
        """ Sum counts across all samples for kraken data """

        # Sum the percentages for each taxa across all samples
        # Allows us to pick top-5 for each rank
        # Use percentages instead of counts so that deeply-sequences samples
        # are not unfairly over-represented
        for s_name, data in self.kraken_raw_data.items():
            total_guess_count = None
            for row in data:

                # Convenience vars that are easier to read
                rank_code = row['rank_code']
                classif = row['classif']

                # Skip anything that doesn't exactly fit a tax rank level
                if row['rank_code'] == '-' or any(c.isdigit() for c in row['rank_code']):
                    continue

                # Calculate the total read count using percentages
                # We use either unclassified or the first domain encountered, to try to use the largest proportion of reads = most accurate guess
                if rank_code == 'U' or (rank_code == 'D' and row['counts_rooted'] > total_guess_count):
                    self.kraken_sample_total_readcounts[s_name] = round(float(row['counts_rooted']) / (row['percent'] / 100.0))
                    total_guess_count = row['counts_rooted']

                if rank_code not in self.kraken_total_pct:
                    self.kraken_total_pct[rank_code] = dict()
                    self.kraken_total_counts[rank_code] = dict()

                if classif not in self.kraken_total_pct[rank_code]:
                    self.kraken_total_pct[rank_code][classif] = 0
                    self.kraken_total_counts[rank_code][classif] = 0
                self.kraken_total_pct[rank_code][classif] += row['percent']
                self.kraken_total_counts[rank_code][classif] += row['counts_rooted']

    def top_five_barplot(self):
        """ Add a bar plot showing the top-5 from each taxa rank """

        t_ranks = OrderedDict()
        t_ranks['S'] = 'Species'
        t_ranks['G'] = 'Genus'
        t_ranks['F'] = 'Family'
        t_ranks['O'] = 'Order'
        t_ranks['C'] = 'Class'
        t_ranks['P'] = 'Phylum'
        t_ranks['K'] = 'Kingdom'
        t_ranks['D'] = 'Domain'
        t_ranks['R'] = 'Root'
        # t_ranks['U'] = 'Unclassified'

        pd = []
        cats = list()
        pconfig = {
            'id': 'kraken-topfive-plot',
            'title': 'Kraken 2: Top taxa',
            'ylab': 'Number of fragments',
            'data_labels': list(t_ranks.values())
        }

        for rank_code, rank_name in t_ranks.items():
            rank_cats = OrderedDict()
            rank_data = dict()

            # Loop through the summed tax percentages to get the top 5 across all samples
            try:
                sorted_pct = sorted(self.kraken_total_pct[rank_code].items(), key=lambda x: x[1], reverse=True)
            except KeyError:
                # Taxa rank not found in this sample
                continue
            i = 0
            counts_shown = {}
            for classif, pct_sum in sorted_pct:
                i += 1
                if i > 5:
                    break
                rank_cats[classif] = {'name': classif}
                # Pull out counts for this rank + classif from each sample
                for s_name, d in self.kraken_raw_data.items():
                    if s_name not in rank_data:
                        rank_data[s_name] = dict()
                    if s_name not in counts_shown:
                        counts_shown[s_name] = 0
                    for row in d:
                        if row['rank_code'] == rank_code:
                            if row['classif'] == classif:
                                if classif not in rank_data[s_name]:
                                    rank_data[s_name][classif] = 0
                                rank_data[s_name][classif] += row['counts_rooted']
                                counts_shown[s_name] += row['counts_rooted']

            # Add in unclassified reads and "other" - we presume from other species etc.
            for s_name, d in self.kraken_raw_data.items():
                for row in d:
                    if row['rank_code'] == 'U':
                        rank_data[s_name]['U'] = row['counts_rooted']
                        counts_shown[s_name] += row['counts_rooted']
                rank_data[s_name]['other'] = self.kraken_sample_total_readcounts[s_name] - counts_shown[s_name]

                # This should never happen... But it does sometimes if the total read count is a bit off
                if rank_data[s_name]['other'] < 0:
                    log.debug("Found negative 'other' count for {} ({}): {}".format(s_name, t_ranks[rank_code], rank_data[s_name]['other']))
                    rank_data[s_name]['other'] = 0

            rank_cats['other'] = { 'name': 'Other', 'color': '#cccccc' }
            rank_cats['U'] = { 'name': 'Unclassified', 'color': '#d4949c' }

            cats.append(rank_cats)
            pd.append(rank_data)

        self.add_section (
            name = 'Top taxa',
            anchor = 'kraken-topfive',
            description = 'The number of reads falling into the top 5 taxa across different ranks.',
            helptext = """
                To make this plot, the percentage of each sample assigned to a given taxa is summed across all samples.
                The counts for these top five taxa are then plotted for each of the 9 different taxa ranks.
                The unclassified count is always shown across all taxa ranks.

                The total number of reads is approximated by dividing the number of `unclassified` reads by the percentage of
                the library that they account for.
                Note that this is only an approximation, and that kraken percentages don't always add to exactly 100%.

                The category _"Other"_ shows the difference between the above total read count and the sum of the read counts
                in the top 5 taxa shown + unclassified. This should cover all taxa _not_ in the top 5, +/- any rounding errors.

                Note that any taxon that does not exactly fit a taxon rank (eg. `-` or `G2`) is ignored.
            """,
            plot = bargraph.plot(pd, cats, pconfig)
        )
