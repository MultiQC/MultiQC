#!/usr/bin/env python

""" MultiQC module to parse output from HiCUP """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='HiCUP', anchor='hicup',
        href='http://www.bioinformatics.babraham.ac.uk/projects/hicup/', 
        info="(Hi-C User Pipeline) is tool for mapping and performing "\
         "quality control on Hi-C data.")
        
        # Find and load any HiCUP summary reports
        self.hicup_data = dict()
        for f in self.find_log_files(config.sp['hicup']):
            self.parse_hicup_logs(f)

        if len(self.hicup_data) == 0:
            log.debug("Could not find any HiCUP data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.hicup_data)))
        
        import json
        print(json.dumps(self.hicup_data, indent=4))
        
        # Write parsed data to a file
        self.write_data_file(self.hicup_data, 'multiqc_hicup')
        
        # Basic Stats Table
        self.hicup_stats_table()
        
        # Report sections
        self.sections = list()
        # self.sections.append({
        #     'name': 'Truncating and Mapping',
        #     'anchor': 'hicup-truncating-mapping',
        #     'content': self.hicup_truncating_chart() + self.hicup_alignment_chart()
        # })
        # 
        # self.sections.append({
        #     'name': 'Filtering',
        #     'anchor': 'hicup-filtering',
        #     'content': self.hicup_filtering_chart()
        # })
        # 
        # self.sections.append({
        #     'name': 'Di-tag length Distribution',
        #     'anchor': 'hicup-lengths',
        #     'content': self.hicup_lengths_chart()
        # })
        # 
        # self.sections.append({
        #     'name': 'De-Duplication',
        #     'anchor': 'hicup-deduplication',
        #     'content': self.hicup_dedup_chart()
        # })


    def parse_hicup_logs(self, f):
        """ Parse a HiCUP summary report """
        if not f['fn'].endswith('.txt'):
            return None
        header = []
        lines = f['f'].splitlines()
        for l in lines:
            s = l.split("\t")
            if len(header) == 0:
                if s[0] != 'File':
                    return None
                header = s[1:]
            else:
                s_name = s[0].lstrip('HiCUP_output/')
                parsed_data = {}
                for idx, num in enumerate(s[1:]):
                    try:
                        parsed_data[header[idx]] = float(num)
                    except:
                        parsed_data[header[idx]] = num
                if s_name in self.hicup_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.add_data_source(f, s_name)
                self.hicup_data[s_name] = parsed_data


    def hicup_stats_table(self):
        """ Add core HiCUP stats to the general stats table """
        headers = OrderedDict()
        headers['Percentage_Ditags_Passed_Through_HiCUP'] = {
            'title': '% Passed',
            'description': 'Percentage Di-Tags Passed Through HiCUP',
            'max': 100,
            'min': 0,
            'scale': 'YlGn',
            'format': '{:.1f}%',
        }
        headers['Deduplication_Read_Pairs_Uniques'] = {
            'title': 'M Unique',
            'description': 'Unique Di-Tags (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['Percentage_Uniques'] = {
            'title': '% Duplicates',
            'description': 'Percent Duplicate Di-Tags',
            'max': 100,
            'min': 0,
            'scale': 'YlGn-rev',
            'modify': lambda x: 100 - x,
            'format': '{:.1f}%',
        }
        headers['Valid_Pairs'] = {
            'title': 'M Valid',
            'description': 'Valid Pairs (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['Percentage_Valid'] = {
            'title': '% Valid',
            'description': 'Percent Valid Pairs',
            'max': 100,
            'min': 0,
            'scale': 'YlGn',
            'format': '{:.1f}%',
        }
        headers['Paired_Read_1'] = {
            'title': 'M Pairs Aligned',
            'description': 'Paired Alignments (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['Percentage_Mapped'] = {
            'title': '% Aligned',
            'description': 'Percentage of Paired Alignments',
            'max': 100,
            'min': 0,
            'scale': 'YlGn',
            'format': '{:.1f}%',
        }
        self.general_stats_addcols(self.hicup_data, headers, 'hicup')
    
    def hicup_length_trimmed_plot (self):
        """ Generate the hicup plot """    
        pconfig = {
            'id': 'hicup_plot',
            'title': 'hicup complexity curve',
            'ylab': 'Unique Molecules',
            'xlab': 'Total Molecules (including duplicates)',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} total</b>: {point.y:,.0f} unique',
            'extra_series': [{
                'name': 'x = y',
                'data': [[0, 0], [self.total_max, self.total_max]],
                'dashStyle': 'Dash',
                'lineWidth': 1,
                'color': '#000000',
                'marker': { 'enabled': False },
                'enableMouseTracking': False,
                'showInLegend': False,
            }]
        }
        return "<p>A shallow curve indicates complexity saturation. The dashed line \
                shows a perfectly complex library where total reads = unique reads.</o>" \
                 + self.plot_xy_data(self.hicup_data, pconfig)
