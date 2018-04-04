#!/usr/bin/env python

""" MultiQC module to parse output from HiCPro """

from __future__ import print_function
from collections import OrderedDict
import os.path
import logging

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ HiCPro module, parses log and stats files saved by HiCPro. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='HiCPro', anchor='hicpro',
        href='https://github.com/nservant/HiC-Pro',
        info="HiC-Pro is an efficient and flexible pipeline for Hi-C data processing")

        log.info("Hello World!")
        
        # Find and load any HiCPro summary reports
        self.hicpro_data = dict()
        for f in self.find_log_files('hicpro'):
            self.parse_hicpro_stats(f)

        # Filter to strip out ignored sample names
        self.hicpro_data = self.ignore_samples(self.hicpro_data)

        if len(self.hicpro_data) == 0:
            raise UserWarning

        # Find all files for mymod
        #for myfile in self.find_log_files('hicpro'):
        #    print( myfile['fn'] )      # Filename    
        #    print( myfile['s_name'] )  # Sample name (from cleaned filename)
        #    print( myfile['root'] )    # Directory file was in
    
        log.info("Found {} reports".format(len(self.hicpro_data)))
        
        # Write parsed data to a file
        #self.write_data_file(self.hicup_data, 'multiqc_hicup')

        # Basic Stats Table
        self.hicpro_stats_table()

        # Report sections
        self.add_section (
            name = 'Read Mapping',
            anchor = 'hicpro-mapping',
            plot = self.hicpro_mapping_chart()
        )

        self.add_section (
            name = 'Read Pair Filtering',
            anchor = 'hicup-filtering',
            plot = self.hicpro_filtering_chart()
        )

        self.add_section (
            name = 'Duplicates removal',
            anchor = 'hicpro-rmdup',
            plot = self.hicpro_rmdup_chart()
        )

        self.add_section (
            name = 'Allele-specific analysis',
            anchor = 'hicpro-asan',
            plot = self.hicpro_as_chart()
        )
        
        self.add_section (
             name = 'Length Distribution',
             anchor = 'hicpro-lengths',
             plot = self.hicpro_insertsize_chart()
        )

        

    def parse_hicpro_stats(self, f):
        """ Parse a HiCPro stat file """
        s_name = os.path.basename(f['root'])
        data = {}
        for l in f['f'].splitlines():
            if not l.startswith('#'):
                s = l.split("\t")
                data[s[0]] = s[1]
                
        if s_name in self.hicpro_data:
            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))

        self.add_data_source(f, s_name)
        self.hicpro_data[s_name] = data


    def hicpro_stats_table(self):
         """ Add HiC-Pro stats to the general stats table """
         headers = OrderedDict()
    #     headers['Percentage_Ditags_Passed_Through_HiCUP'] = {
    #         'title': '% Passed',
    #         'description': 'Percentage Di-Tags Passed Through HiCUP',
    #         'max': 100,
    #         'min': 0,
    #         'suffix': '%',
    #         'scale': 'YlGn'
    #     }
    #     headers['Deduplication_Read_Pairs_Uniques'] = {
    #         'title': '{} Unique'.format(config.read_count_prefix),
    #         'description': 'Unique Di-Tags ({})'.format(config.read_count_desc),
    #         'min': 0,
    #         'scale': 'PuRd',
    #         'modify': lambda x: x * config.read_count_multiplier,
    #         'shared_key': 'read_count'
    #     }
    #     headers['Percentage_Uniques'] = {
    #         'title': '% Duplicates',
    #         'description': 'Percent Duplicate Di-Tags',
    #         'max': 100,
    #         'min': 0,
    #         'suffix': '%',
    #         'scale': 'YlGn-rev',
    #         'modify': lambda x: 100 - x
    #     }
         headers['Valid_Pairs'] = {
             'title': '{} Valid'.format(config.read_count_prefix),
             'description': 'Valid Pairs ({})'.format(config.read_count_desc),
             'min': 0,
             'scale': 'PuRd',
             'modify': lambda x: x * config.read_count_multiplier,
             'shared_key': 'read_count'
         }
    #     headers['Percentage_Valid'] = {
    #         'title': '% Valid',
    #         'description': 'Percent Valid Pairs',
    #         'max': 100,
    #         'min': 0,
    #         'suffix': '%',
    #         'scale': 'YlGn'
    #     }
    #     headers['Paired_Read_1'] = {
    #         'title': '{} Pairs Aligned'.format(config.read_count_prefix),
    #         'description': 'Paired Alignments ({})'.format(config.read_count_desc),
    #         'min': 0,
    #         'scale': 'PuRd',
    #         'modify': lambda x: x * config.read_count_multiplier,
    #         'shared_key': 'read_count'
    #     }
         headers['Percentage_Mapped'] = {
             'title': '% Aligned',
             'description': 'Percentage of Paired Alignments',
             'max': 100,
             'min': 0,
             'suffix': '%',
             'scale': 'YlGn'
         }
         self.general_stats_addcols(self.hicpro_data, headers, 'HiC-Pro')

    def hicpro_mapping_chart (self):
        """ Generate the HiCUP Aligned reads plot """

        # Specify the order of the different possible categories
    #     keys = OrderedDict()
    #     keys['Unique_Alignments_Read']   = { 'color': '#2f7ed8', 'name': 'Unique Alignments' }
    #     keys['Multiple_Alignments_Read'] = { 'color': '#492970', 'name': 'Multiple Alignments' }
    #     keys['Failed_To_Align_Read']     = { 'color': '#0d233a', 'name': 'Failed To Align' }
    #     keys['Too_Short_To_Map_Read']    = { 'color': '#f28f43', 'name': 'Too short to map' }

    #     # Construct a data structure for the plot - duplicate the samples for read 1 and read 2
    #     data = {}
    #     for s_name in self.hicup_data:
    #         data['{} Read 1'.format(s_name)] = {}
    #         data['{} Read 2'.format(s_name)] = {}
    #         data['{} Read 1'.format(s_name)]['Unique_Alignments_Read'] = self.hicup_data[s_name]['Unique_Alignments_Read_1']
    #         data['{} Read 2'.format(s_name)]['Unique_Alignments_Read'] = self.hicup_data[s_name]['Unique_Alignments_Read_2']
    #         data['{} Read 1'.format(s_name)]['Multiple_Alignments_Read'] = self.hicup_data[s_name]['Multiple_Alignments_Read_1']
    #         data['{} Read 2'.format(s_name)]['Multiple_Alignments_Read'] = self.hicup_data[s_name]['Multiple_Alignments_Read_2']
    #         data['{} Read 1'.format(s_name)]['Failed_To_Align_Read'] = self.hicup_data[s_name]['Failed_To_Align_Read_1']
    #         data['{} Read 2'.format(s_name)]['Failed_To_Align_Read'] = self.hicup_data[s_name]['Failed_To_Align_Read_2']
    #         data['{} Read 1'.format(s_name)]['Too_Short_To_Map_Read'] = self.hicup_data[s_name]['Too_Short_To_Map_Read_1']
    #         data['{} Read 2'.format(s_name)]['Too_Short_To_Map_Read'] = self.hicup_data[s_name]['Too_Short_To_Map_Read_2']

    #     # Config for the plot
    #     config = {
    #         'id': 'hicup_mapping_stats_plot',
    #         'title': 'HiCUP: Mapping Statistics',
    #         'ylab': '# Reads',
    #         'cpswitch_counts_label': 'Number of Reads'
    #     }

    #     return bargraph.plot(data, keys, config)

    def hicpro_filtering_chart (self):
        """ Generate the HiC-Pro filtering plot """

    #     # Specify the order of the different possible categories
    #     keys = OrderedDict()
    #     keys['Valid_Pairs'] =            { 'color': '#2f7ed8', 'name': 'Valid Pairs' }
    #     keys['Same_Fragment_Internal'] = { 'color': '#0d233a', 'name': 'Same Fragment - Internal' }
    #     keys['Same_Circularised'] =      { 'color': '#910000', 'name': 'Same Fragment - Circularised' }
    #     keys['Same_Dangling_Ends'] =     { 'color': '#8bbc21', 'name': 'Same Fragment - Dangling Ends' }
    #     keys['Re_Ligation'] =            { 'color': '#1aadce', 'name': 'Re-ligation' }
    #     keys['Contiguous_Sequence'] =    { 'color': '#f28f43', 'name': 'Contiguous Sequence' }
    #     keys['Wrong_Size'] =             { 'color': '#492970', 'name': 'Wrong Size' }

    #     # Config for the plot
    #     config = {
    #         'id': 'hicup_filtering_plot',
    #         'title': 'HiCUP: Filtering Statistics',
    #         'ylab': '# Read Pairs',
    #         'cpswitch_counts_label': 'Number of Read Pairs',
    #         'cpswitch_c_active': False
    #     }

    #     return bargraph.plot(self.hicup_data, keys, config)

    def hicpro_rmdup_chart (self):
        """ Generate the HiC-Pro interaction plot """

    #     # Specify the order of the different possible categories
    #     keys = OrderedDict()
    #     keys['Deduplication_Cis_Close_Uniques'] = { 'color': '#2f7ed8', 'name': 'Unique: cis < 10Kbp' }
    #     keys['Deduplication_Cis_Far_Uniques']   = { 'color': '#0d233a', 'name': 'Unique: cis > 10Kbp' }
    #     keys['Deduplication_Trans_Uniques']     = { 'color': '#492970', 'name': 'Unique: trans' }
    #     keys['Duplicate_Read_Pairs']            = { 'color': '#f28f43', 'name': 'Duplicate read pairs' }

    #     # Config for the plot
    #     config = {
    #         'id': 'hicup_dedup_plot',
    #         'title': 'HiCUP: De-Duplication Statistics',
    #         'ylab': '# Di-Tags',
    #         'cpswitch_counts_label': 'Number of Di-Tags',
    #         'cpswitch_c_active': False
    #     }

    #     return bargraph.plot(self.hicup_data, keys, config)

    def hicpro_as_chart (self):
        """ Generate Allele-specific plot"""

    def hicpro_insertsize_chart (self):
        """ Generate insert size histogram"""



        
