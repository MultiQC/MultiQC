#!/usr/bin/env python

""" MultiQC module to parse output from CRISPResso """

from __future__ import print_function
from collections import OrderedDict
import logging
from pprint import pprint
import os
import re

import pandas as pd

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ CRISPResso module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='CRISPResso2', anchor='crispresso2',
        href="http://crispresso2.pinellolab.org",
        info="Analysis of genome editing outcomes from deep sequencing data.")

        # Find and load any Kallisto reports
        self.crispresso2_data = dict()
        for f in self.find_log_files('crispresso2/allele_frequency_table', filehandles=True):
            # log.info("Found {}!")
            filename = os.path.join(f['root'], f['fn'])
            # pprint(f)
            self.parse_allele_frequencies(filename)

        for f in self.find_log_files('crispresso2/mapping_statistics'):
            # log.info("Found {}!")
            # filename = os.path.join(f['root'], f['fn'])
            pprint(f)
            self.parse_mapping_statistics(f)

        for f in self.find_log_files('crispresso2/editing_frequencies'):
            # log.info("Found {}!")
            # filename = os.path.join(f['root'], f['fn'])
            pprint(f)


        # # Filter to strip out ignored sample names
        # self.kallisto_data = self.ignore_samples(self.kallisto_data)
        #
        # if len(self.kallisto_data) == 0:
        #     raise UserWarning
        #
        log.info("Found {} reports".format(len(self.crispresso2_data)))

        # Write parsed report data to a file (BaseMultiqcModule function)
        self.write_data_file(self.crispresso2_data, 'multiqc_crispresso2')
        pprint(self.crispresso2_data)

        #
        # Plot of alleles
        self.add_section( plot = self.crispresso2_mapping_stats_plot() )
        self.add_section( plot =  self.crispresso2_allele_plot() )
        self.crispresso2_general_stats_table()

    def parse_allele_frequencies(self, filename):
        df = pd.read_csv(filename, sep='\t', compression='zip')
        sample_id = self.get_sample_id_from_dirname(filename)

        allele_data = self.get_read_status_line(df, sample_id, "#Reads")

        # Assign to this sample_id
        self.crispresso2_data[sample_id] = allele_data

    def parse_mapping_statistics(self, f):
        sample_id = self.get_sample_id_from_dirname(f['root'])

        for i, line in enumerate(f['f'].splitlines()):
            split = line.split('\t')
            if i == 0:
                keys = split
            if i == 1:
                values = [int(x) for x in split]

        data = dict(zip(keys, values))
        data['reads_unaligned'] = data['READS AFTER PREPROCESSING'] - data['READS ALIGNED']
        data['reads_removed_by_preprocessing'] = data['READS IN INPUTS'] - data['READS AFTER PREPROCESSING']

        self.crispresso2_data[sample_id].update(data)

    @staticmethod
    def get_sample_id_from_dirname(dirname):
        # sample_id = os.path.basename(os.path.dirname(filename))
        sample_id = dirname.split("CRISPResso_on_")[-1]
        sample_id = sample_id.split('.collapsed')[0]
        return sample_id

    @staticmethod
    def get_read_status_line(df, sample_id, read_quantification,
                             groupby=['Reference_Name', 'Read_Status']):
        reads = df.groupby(groupby)[read_quantification].sum()
        reads_tidy = reads.reset_index()
        reads_tidy['sample_id'] = sample_id
        reads_tidy['reference_name_with_status'] = reads_tidy.apply(
            lambda x: "{Reference_Name}_{Read_Status}".format(**x), axis=1)
        df_line = reads_tidy.pivot(index='sample_id',
                                   columns='reference_name_with_status',
                                   values=read_quantification)
        series = df_line.iloc[0]
        percentages = 100 * series/series.sum()
        percentages.index = percentages.index + "_percentage"

        # Combine both percentages and raw numbers
        line = pd.concat([percentages, series])

        # Convert to dictionary since that's what MultiQC uses
        line = line.to_dict()
        return line

    # def kallisto_general_stats_table(self):
    #     """Take the parsed stats from the Kallisto report and add it to the
    #     basic stats table at the top of the report """
    #
    #     headers = OrderedDict()
    #     headers['fragment_length'] = {
    #         'title': 'Frag Length',
    #         'description': 'Estimated average fragment length',
    #         'min': 0,
    #         'suffix': 'bp',
    #         'scale': 'RdYlGn'
    #     }
    #     headers['percent_aligned'] = {
    #         'title': '% Aligned',
    #         'description': '% processed reads that were pseudoaligned',
    #         'max': 100,
    #         'min': 0,
    #         'suffix': '%',
    #         'scale': 'YlGn'
    #     }
    #     headers['pseudoaligned_reads'] = {
    #         'title': '{} Aligned'.format(config.read_count_prefix),
    #         'description': 'Pseudoaligned reads ({})'.format(config.read_count_desc),
    #         'min': 0,
    #         'scale': 'PuRd',
    #         'modify': lambda x: x * config.read_count_multiplier,
    #         'shared_key': 'read_count'
    #     }
    #     self.general_stats_addcols(self.kallisto_data, headers)

    def crispresso2_allele_plot(self):
        """ Make the HighCharts HTML to plot the alignment rates """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['HDR_MODIFIED'] = { 'color': "#6a3d9a", 'name': 'HDR Modified' }
        keys['HDR_UNMODIFIED'] = { 'color': "#cab2d6", 'name': 'HDR Unmodified' }
        keys['Reference_MODIFIED'] = { 'color': "#33a02c", 'name': 'Reference Modified' }
        keys['Reference_UNMODIFIED'] = { 'color': "#b2df8a", 'name': 'Reference Unmodified' }

        # Config for the plot
        config = {
            'id': 'crispresso2_alleles',
            'title': 'CRISPResso2: Allele Frequencies',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(self.crispresso2_data, keys, config)

    def crispresso2_mapping_stats_plot(self):
        """ Make the HighCharts HTML to plot the alignment rates """

        # Specify the order of the different possible categories

        ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
         "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"]
        keys = OrderedDict()
        keys['reads_removed_by_preprocessing'] = { 'name': 'Removed by preprocessing' }
        keys['READS ALIGNED'] = {  'name': 'Aligned' }
        keys['reads_unaligned'] = { 'name': 'Unaligned' }

        # Config for the plot
        config = {
            'id': 'crispresso2_mapping_stats',
            'title': 'CRISPResso2: Mapping Statistics',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(self.crispresso2_data, keys, config)

    def crispresso2_general_stats_table(self):
        """ Take the parsed stats from the CRISPresso2 report and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()

        # headers['READS IN INPUTS'] = {
        #     'title': 'Input reads to CRISPResso2',
        #     'min': 0,
        #     # 'suffix': '%',
        #     # 'format': '{:,d}',
        #     # 'scale': 'Blues'
        # }
        # headers['READS AFTER PREPROCESSING'] = {
        #     'title': 'Reads after CRISPResso2 preprocessing',
        #     'min': 0,
        #     # 'suffix': '%',
        #     # 'format': '{:,d}',
        #     # 'scale': 'Blues'
        # }
        # headers['READS ALIGNED'] = {
        #     'title': 'Reads aligned in CRISPResso2',
        #     'min': 0,
        #     # 'suffix': '%',
        #     # 'format': '{:,d}',
        #     # 'scale': 'Blues'
        # }
        headers['HDR_UNMODIFIED_percentage'] = {
            'title': 'HDR Unmodified',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:,.2f}',
            'scale': 'Purples'
        }
        headers['HDR_MODIFIED_percentage'] = {
            'title': 'HDR Modified',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:,.2f}',
            'scale': 'Purples'
        }
        headers['Reference_UNMODIFIED_percentage'] = {
            'title': 'Reference Modified',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:,.2f}',
            'scale': 'Greens'
        }
        headers['Reference_MODIFIED_percentage'] = {
            'title': 'Reference Un,odified',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:,.2f}',
            'scale': 'Greens'
        }
        self.general_stats_addcols(self.crispresso2_data, headers)
