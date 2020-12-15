#!/usr/bin/env python

""" MultiQC module to parse output from pbmarkdup"""

import logging

from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    pbmarkdup module class, parses pbmarkdup output.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
                name='pbmarkdup', anchor='pbmarkdup',
                href='https://github.com/PacificBiosciences/pbmarkdup',
                info=("""
                 takes one or multiple sequencing chips of an amplified libray
                 as HiFi reads and marks or removes duplicates.
                 """)
        )

        self.pbmarkdup = dict()

        for logfile in self.find_log_files('pbmarkdup', filehandles=True):
            self.pbmarkdup[logfile['s_name']] = self.parse_logfile(logfile)

        # Filter to strip out ignored sample names
        self.pbmarkdup = self.ignore_samples(self.pbmarkdup)

        # Raise UserWarning if we did not find any data
        if not self.pbmarkdup:
            raise UserWarning

        # Write the parsed data to file
        self.write_data_file(self.pbmarkdup, 'multiqc_pbmarkdup')

        self.pbmarkdup_add_general_stats()

    def parse_logfile(self, logfile):
        """
        Parse the standard output from pbmarkdup

        This is currently quite ugly, since pbmarkdup does not have an easily
        parsable output format. See this github issue:
        https://github.com/PacificBiosciences/pbbioconda/issues/365

        For now, we just assume that the data fields are always in the same
        order, without any missing values.
        """

        file_content = logfile['f']

        # This is not tab separated :(
        expected_header = 'LIBRARY         READS    UNIQUE MOLECULES    DUPLICATE READS'
        header = next(file_content).strip()

        assert header == expected_header, f'Unknown header: "{header}"'

        data = dict()

        # Each parsable line is either for a library, or for 'TOTAL', the sum
        # of all parsed libraries
        for line in file_content:
            # Lines which start with --- denote separators in the file, and do
            # not need to be parsed
            if line.startswith('-'):
                continue

            # Not very nice, we assume that all fields are always present
            name, reads, unique_mol_count, unique_mol_perc, duplicate_count, duplicate_perc = line.split()

            # We are only interested in the counts, not the percentages
            data[name] = {
                    'READS': int(reads),
                    'UNIQUE MOLECULES': int(unique_mol_count),
                    'DUPLICATE READS': int(duplicate_count)
            }

        return data

    def pbmarkdup_add_general_stats(self):
        """ Add pbmarkdup duplicates to the general stats table """

        general_stats_headers = OrderedDict([
            ('unique_molecules', {
                'id': 'unique_molecules',
                'title': 'Unique Molecules',
                'description': 'Percentage of unique molecules',
                'suffix': '%',
                'min': 0,
                'max': 100,
                'modify': lambda x: x * 100,
                'scale': 'RdYlGn'
                }
            ),
            ('duplicate_reads', {
                'id': 'duplicate_erads',
                'title': 'Duplicate Reads',
                'description': 'Percentage of duplicate reads',
                'suffix': '%',
                'min': 0,
                'max': 100,
                'modify': lambda x: x * 100,
                'scale': 'RdYlGn-rev'
                }
            )
        ])

        general = dict()

        for sample, sample_stats in self.pbmarkdup.items():
            # We are only interested in the total counts per sample
            total = sample_stats['TOTAL']

            # Calculate the percentages unique and duplicated
            perc_unique = total['UNIQUE MOLECULES']/total['READS']
            perc_duplicate = total['DUPLICATE READS']/total['READS']

            general[sample] = {
                    'unique_molecules': perc_unique,
                    'duplicate_reads': perc_duplicate
            }

        self.general_stats_addcols(general, general_stats_headers)
