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
        for myfile in self.find_log_files('hicpro'):
            print( myfile['fn'] )      # Filename    
            print( myfile['s_name'] )  # Sample name (from cleaned filename)
            print( myfile['root'] )    # Directory file was in
    
        log.info("Found {} reports".format(len(self.hicpro_data)))
        
        # Write parsed data to a file
        #self.write_data_file(self.hicup_data, 'multiqc_hicup')

        # Basic Stats Table
        self.hicpro_stats_table()

        # Report sections
        self.add_section (
            name = 'Read Mapping',
            anchor = 'hicpro-mapping',
            description = 'todo - description',
            helptext = 'todo - help',
            plot = self.hicpro_mapping_chart()
        )

        self.add_section (
            name = 'Read Pairing',
            anchor = 'hicpro-pairing',
            plot = self.hicpro_pairing_chart()
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

        ## Check if data already exists
        if s_name in self.hicpro_data.keys():
            data = self.hicpro_data[s_name]
        else:
            data = {}

        for l in f['f'].splitlines():
            if not l.startswith('#'):
                s = l.split("\t")
                data[s[0]] = int(s[1])
                
        if s_name in self.hicpro_data:
            log.debug("Duplicated sample name found! Overwriting: {}".format(s_name))

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
         headers['Valid_interaction_pairs'] = {
             'title': '{} Valid', #.format(config.read_count_prefix),
             'description': 'Valid Pairs ({})', #.format(config.read_count_desc),
             'min': 0,
             'scale': 'PuRd',
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
        """ Generate the HiC-Pro Aligned reads plot """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['Full_Alignments_Read']   = { 'color': '#2f7ed8', 'name': 'Full reads Alignments' }
        keys['Trimmed_Alignments_Read'] = { 'color': '#492970', 'name': 'Trimmed reads Alignments' }
        keys['Failed_To_Align_Read']     = { 'color': '#0d233a', 'name': 'Failed To Align' }
    
        # Construct a data structure for the plot - duplicate the samples for read 1 and read 2
        data = {}
        for s_name in self.hicpro_data:
            data['{} Read 1'.format(s_name)] = {}
            data['{} Read 2'.format(s_name)] = {}
            data['{} Read 1'.format(s_name)]['Full_Alignments_Read'] = self.hicpro_data[s_name]['global']
            data['{} Read 2'.format(s_name)]['Full_Alignments_Read'] = self.hicpro_data[s_name]['global']
            data['{} Read 1'.format(s_name)]['Trimmed_Alignments_Read'] = self.hicpro_data[s_name]['local']
            data['{} Read 2'.format(s_name)]['Trimmed_Alignments_Read'] = self.hicpro_data[s_name]['local']
            data['{} Read 1'.format(s_name)]['Failed_To_Align_Read'] = int(self.hicpro_data[s_name]['total']) - int(self.hicpro_data[s_name]['mapped'])
            data['{} Read 2'.format(s_name)]['Failed_To_Align_Read'] = int(self.hicpro_data[s_name]['total']) - int(self.hicpro_data[s_name]['mapped'])

        # Config for the plot
        config = {
            'id': 'hicpro_mapping_stats_plot',
            'title': 'HiC-Pro: Mapping Statistics',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(data, keys, config)

    def hicpro_pairing_chart (self):
        """ Generate Pairing chart """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['Unique_Pairs'] = { 'color': '#2f7ed8', 'name': 'Uniquely Aligned' }
        keys['Low_Qual_Pairs'] = { 'color': '#492970', 'name': 'Low Quality' }
        keys['Singleton_Pairs'] = { 'color': '#0d233a', 'name': 'Singleton' }
        keys['Multi_Pairs'] = { 'color': '#0d233a', 'name': 'Multi Aligned' }
        keys['Failed_To_Align_Pairs'] = { 'color': '#0d233a', 'name': 'Failed To Align' }
    
        # Construct a data structure for the plot - duplicate the samples for read 1 and read 2
        data = {}
        for s_name in self.hicpro_data:
            print ( self.hicpro_data[s_name] )
            data['{} Pairs'.format(s_name)] = {}
            data['{} Pairs'.format(s_name)]['Failed_To_Align_Pairs'] = self.hicpro_data[s_name]['Unmapped_pairs']
            data['{} Pairs'.format(s_name)]['Low_Qual_Pairs'] = self.hicpro_data[s_name]['Low_qual_pairs']
            data['{} Pairs'.format(s_name)]['Unique_Pairs'] = self.hicpro_data[s_name]['Unique_paired_alignments']
            data['{} Pairs'.format(s_name)]['Multi_Pairs'] = self.hicpro_data[s_name]['Multiple_pairs_alignments']
            data['{} Pairs'.format(s_name)]['Singleton_Pairs'] = self.hicpro_data[s_name]['Pairs_with_singleton']

        # Config for the plot
        config = {
        'id': 'hicpro_pairing_stats_plot',
            'title': 'HiC-Pro: Pairing Statistics',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(data, keys, config)

        
    def hicpro_filtering_chart (self):
        """ Generate the HiC-Pro filtering plot """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['Valid_interaction_pairs_FF'] = { 'name': 'Valid Pairs FF' }
        keys['Valid_interaction_pairs_RR'] = { 'name': 'Valid Pairs RR' }
        keys['Valid_interaction_pairs_RF'] = { 'name': 'Valid Pairs RF' }
        keys['Valid_interaction_pairs_FR'] = { 'name': 'Valid Pairs FR' }
        keys['Self_Circle_pairs'] = { 'name': 'Same Fragment - Self-Circle' }
        keys['Dangling_ends_pairs'] = { 'name': 'Same Fragment - Dangling Ends' }
        keys['Religation_pairs'] = { 'name': 'Re-ligation' }
        keys['Dumped_pairs'] = { 'name': 'Dumped pairs' }

        # Config for the plot
        config = {
            'id': 'hicpro_filtering_plot',
            'title': 'HiC-Pro: Filtering Statistics',
            'ylab': '# Read Pairs',
            'cpswitch_counts_label': 'Number of Read Pairs',
            'cpswitch_c_active': False
        }
        
        return bargraph.plot(self.hicpro_data, keys, config)

    def hicpro_rmdup_chart (self):
        """ Generate the HiC-Pro interaction plot """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['cis_shortRange'] = { 'color': '#2f7ed8', 'name': 'Unique: cis < 10Kbp' }
        keys['cis_longRange'] = { 'color': '#0d233a', 'name': 'Unique: cis > 10Kbp' }
        keys['trans_interaction']  = { 'color': '#492970', 'name': 'Unique: trans' }
        keys['valid_interaction_rmdup'] = { 'color': '#f28f43', 'name': 'Duplicate read pairs' }

        # Config for the plot
        config = {
            'id': 'hicpro_rmdup_plot',
            'title': 'HiC-Pro: Contact Statistics',
            'ylab': '# Pairs',
            'cpswitch_counts_label': 'Number of Pairs',
            'cpswitch_c_active': False
        }
        
        return bargraph.plot(self.hicpro_data, keys, config)
  
    def hicpro_as_chart (self):
        """ Generate Allele-specific plot"""

    def hicpro_insertsize_chart (self):
        """ Generate insert size histogram"""



        
