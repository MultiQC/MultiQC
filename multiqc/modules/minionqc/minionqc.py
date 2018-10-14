#!/usr/bin/env python

""" MultiQC submodule to parse output from MinIONQC summary stats """

import yaml
import os
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='MinIONQC', anchor='mymod',
        href="https://github.com/roblanf/minion_qc",
        info=" is a QC tool for Oxford Nanopore sequencing data")

        # Find and load any minionqc reports
        self.minionqc_data = dict()
        for f in self.find_log_files('minionqc', filehandles=True):            
            # get sample name
            s_name = self.clean_s_name(os.path.basename(f['root']), f['root'])
            if s_name in self.minionqc_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            
            # adds files used in MultiQC report
            self.add_data_source(f, s_name)

            # parses minionqc summary data
            self.minionqc_data[s_name] = self.parse_minionqc_report(s_name, f['f'])

        # Filter to strip out ignored sample names
        self.minionqc_data = self.ignore_samples(self.minionqc_data)
        if len(self.minionqc_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.minionqc_data)))        

        # columns to present in MultiQC table
        headers = self.headers_to_use()
        self.general_stats_addcols(self.minionqc_data, headers)
        
        # writes parsed data into a file (by MultiQC)
        self.write_data_file(self.minionqc_data, 'multiqc_minionqc')
            
        # plots in order
        self.base_count_plot()
        self.median_q_plot()
        self.reads_n50_plot()

    
    def parse_minionqc_report(self, s_name, f):
        '''
        Parses minionqc's 'summary.yaml' report file for results.
        Uses only the "All reads" stats. Ignores "Q>=x" part.
        '''
        summary_dict = yaml.safe_load(f)

        data_dict = {}
        data_dict = summary_dict['All reads']
        
        # flatten the nested dictionary part 
        for key_1 in ['reads', 'gigabases']:
            for key_2 in data_dict[key_1]:
                data_dict[f'{key_1} {key_2}'] = data_dict[key_1][key_2]
            data_dict.pop(key_1)    # removes key after flattening

        return data_dict


    def headers_to_use(self):
        '''
        Defines features of columns to be used in multiqc table
        '''
        headers = OrderedDict()

        headers['total.reads'] = {
            'title': 'Total reads',
            'description': 'Total number of reads',
            'format': '{:,.0f}',
        }
        headers['total.gigabases'] = {
            'title': 'Total bases (GB)',
            'description': 'Total bases',
            'format': '{:,.2f}',
            # 'suffix': ' GB'
        }
        headers['N50.length'] = {
            'title': 'Reads N50',
            'description': 'Minimum read length needed to cover 50% of all reads',
            'format': '{:,.0f}',
        }
        headers['mean.length'] = {
            'title': 'Mean length',
            'description': 'Mean read length',
            'format': '{:,.0f}',
            'hidden': True,
        }
        headers['median.length'] = {
            'title': 'Median length',
            'description': 'Median read length',
            'format': '{:,.0f}',
        }
        headers['mean.q'] = {
            'title': 'Mean q',
            'description': 'Mean quality of reads',
            'min': 0,
            'max': 20,
            'format': '{:,.1f}',
            'hidden': True,
        }
        headers['median.q'] = {
            'title': 'Median q',
            'description': 'Median quality of reads',
            'min': 0,
            'max': 20,
            'format': '{:,.1f}',
        }
        for s in ['>10kb', '>50kb', '>100kb']:
            headers[f'reads {s}'] = {
                'title': f'Total reads {s}',
                'description': f'Total number of reads {s}',
                'format': '{:,.0f}',
            }
        for s in ['>10kb', '>50kb', '>100kb']:
            headers[f'gigabases {s}'] = {
                'title': f'Bases {s}',
                'description': f'Total bases from reads {s}',
                'format': '{:,.3f}',
                'suffix': ' GB',
                'hidden': True,
            }

        return headers


    def reads_n50_plot (self):
        """ Stacked bar plot showing N50 of reads """
        pconfig = {
            'id': 'minionqc_reads_N50_plot',
            'title': 'MinIONQC: Reads N50',
            'ylab': 'Reads N50',
            'cpswitch': False,
            'yDecimals': False,
            'tt_percentages': False
        }

        # data for plotting
        data = {}
        for sample in self.minionqc_data:
            data[sample] = {}
            data[sample]['N50'] = self.minionqc_data[sample]['N50.length']

        self.add_section (
            name = 'Reads N50',
            anchor = 'minionqc_reads_N50',
            description = 'N50 of reads for each sample.',
            helptext = '''
            Plots N50 of each sample. N50 is the minimum read length needed to cover 50% of all reads.
            ''',

            plot = bargraph.plot(data, pconfig=pconfig)
        )

        return None


    def median_q_plot (self):
        """ Stacked bar plot showing N50 of reads """
        pconfig = {
            'id': 'minionqc_median_q_plot',
            'title': 'MinIONQC: Median q',
            'ylab': 'Median q',
            'cpswitch': False,
            'yDecimals': False,
            'tt_decimals': 1,
            'tt_percentages': False
        }

        # data for plotting
        data = {}
        for sample in self.minionqc_data:
            data[sample] = {}
            data[sample]['Median q'] = self.minionqc_data[sample]['median.q']
        
        self.add_section (
            name = 'Median q',
            anchor = 'minionqc_median_q',
            description = 'Median q of reads for each sample.',
            helptext = '''
            Median read quality for each sample
            ''',

            plot = bargraph.plot(data, pconfig=pconfig)
        )

        return None


    def base_count_plot (self):
        """ Stacked bar plot showing N50 of reads """
        pconfig = {
            'id': 'minionqc_basecount_plot',
            'title': 'MinIONQC: Total base count',
            'ylab': 'Base count (GB)',
            'cpswitch': False,
            'yDecimals': False,
            'tt_decimals': 3,
            'tt_percentages': False
        }

        # data for plotting
        data = {}
        for sample in self.minionqc_data:
            data[sample] = {}
            data[sample]['Base count (GB)'] = self.minionqc_data[sample]['total.gigabases']
        
        self.add_section (
            name = 'Base count',
            anchor = 'minionqc_basecount',
            description = 'Base count for each sample.',
            helptext = 'Total number of bases.',

            plot = bargraph.plot(data, pconfig=pconfig)
        )

        return None
