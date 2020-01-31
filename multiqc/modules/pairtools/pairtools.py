#!/usr/bin/env python

""" MultiQC module to parse stats output from pairtools """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """This MultiQC module parses various
    stats produced by pairtools."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='pairtools', anchor='pairtools',
        href="https://github.com/mirnylab/pairtools",
        info="pairtools is a command-line framework to process sequencing data from a Hi-C experiment.")


        # Find and load any pairtools stats summaries
        self.pairtools_stats = dict()
        for f in self.find_log_files('pairtools', filehandles=True):
            s_name = f['s_name']
            self.pairtools_stats[s_name] = self.parse_pairtools_stats(f)


        # Filter to strip out ignored sample names
        self.pairtools_stats = self.ignore_samples(self.pairtools_stats)

        if len(self.pairtools_stats) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.pairtools_stats)))
        # self.write_data_file(self.pairtools_stats, 'multiqc_pairtools')

        self.pairtools_general_stats()

        # Report sections
        self.add_section (
            name = 'Read Truncation',
            anchor = 'hicup-truncating',
            plot = self.hicup_truncating_chart()
        )



    def parse_pairtools_stats(self, f):
        """ Parse a pairtools summary stats """
        #
        # categories that we need to read from stats file
        # ['chrom_freq',
        #  'cis',
        #  'cis_10kb+',
        #  'cis_1kb+',
        #  'cis_20kb+',
        #  'cis_2kb+',
        #  'cis_40kb+',
        #  'cis_4kb+',
        #  'dist_freq',
        #  'pair_types',
        #  'total',
        #  'total_dups',
        #  'total_mapped',
        #  'total_nodups',
        #  'total_single_sided_mapped',
        #  'total_unmapped',
        #  'trans']
        #
        #
        # s_name = f['s_name']
        # f_name = f['fn']
        # log.info("parsing {} {} ...".format(s_name,f_name))
        #
        # just testing displying random fields from stats ...
        f_handle = f['f']
        _data = dict()
        _i = 0
        num_fields = 8
        for line in f_handle:
            _1,_2 = line.rstrip().split('\t')
            _data[_1] = int(_2)
            _i += 1
            if _i >=num_fields:
                break
        return _data


    def pairtools_general_stats(self):
        """ Add columns to General Statistics table """
        # headers = OrderedDict()
        # headers['total'] = {
        #     'title': 'total',
        #     'description': 'total number of pairs per sample',
        #     'min': 0,
        # }
        # self.general_stats_addcols(self.pairtools_stats, headers, 'pairtools')
        self.general_stats_addcols(self.pairtools_stats)


    def hicup_truncating_chart (self):
        """ Generate the HiCUP Truncated reads plot """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['Not_Truncated_Reads'] = { 'color': '#2f7ed8', 'name': 'Not Truncated' }
        keys['Truncated_Read']      = { 'color': '#0d233a', 'name': 'Truncated' }

        # Construct a data structure for the plot - duplicate the samples for read 1 and read 2
        data = {}
        for s_name in self.hicup_data:
            data['{} Read 1'.format(s_name)] = {}
            data['{} Read 2'.format(s_name)] = {}
            data['{} Read 1'.format(s_name)]['Not_Truncated_Reads'] = self.hicup_data[s_name]['Not_Truncated_Reads_1']
            data['{} Read 2'.format(s_name)]['Not_Truncated_Reads'] = self.hicup_data[s_name]['Not_Truncated_Reads_2']
            data['{} Read 1'.format(s_name)]['Truncated_Read'] = self.hicup_data[s_name]['Truncated_Read_1']
            data['{} Read 2'.format(s_name)]['Truncated_Read'] = self.hicup_data[s_name]['Truncated_Read_2']

        # Config for the plot
        config = {
            'id': 'hicup_truncated_reads_plot',
            'title': 'HiCUP: Truncated Reads',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(data, keys, config)
