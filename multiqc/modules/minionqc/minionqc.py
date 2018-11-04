#!/usr/bin/env python

""" MultiQC submodule to parse output from MinIONQC summary stats """

import yaml
import os
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='MinIONQC', anchor='minionqc',
        href="https://github.com/roblanf/minion_qc",
        info=" is a QC tool for Oxford Nanopore sequencing data")

        # Find and load any minionqc reports
        self.minionqc_data = dict()     # main dataset. Stats from all reads
        self.qfilt_data = dict()        # Stats from quality filtered reads
        q_threshold_list = set()        # quality thresholds
        for f in self.find_log_files('minionqc', filehandles=True):            
            # get sample name
            s_name = self.clean_s_name(os.path.basename(f['root']), f['root'])
            if s_name in self.minionqc_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            
            # adds files used in MultiQC report
            self.add_data_source(f, s_name)

            # parses minionqc summary data
            parsed_dict = self.parse_minionqc_report(s_name, f['f'])
            self.minionqc_data[s_name] = parsed_dict[0]     # stats for all reads
            self.qfilt_data[s_name] = parsed_dict[1]        # stats for q-filtered reads
            q_threshold_list.add(parsed_dict[2])

        # Filter to strip out ignored sample names
        self.minionqc_data = self.ignore_samples(self.minionqc_data)
        if len(self.minionqc_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.minionqc_data)))        

        # columns to present in MultiQC summary table
        headers = self.headers_to_use()
        headers_subset = OrderedDict()
        for k in ['total.gigabases', 'N50.length', 'median.q', 'mean.q']:
            headers_subset[k] = headers[k]
        
        self.general_stats_addcols(self.minionqc_data, headers_subset)
        
        # writes parsed data into a file (by MultiQC)
        self.write_data_file(self.minionqc_data, 'multiqc_minionqc')
            
        # tables and plots in order
        self.table_qALL()
        self.table_qfiltered(q_threshold_list)
        self.base_count_plot()
        self.median_q_plot()
        self.reads_n50_plot()

    
    def parse_minionqc_report(self, s_name, f):
        '''
        Parses minionqc's 'summary.yaml' report file for results.
        Uses only the "All reads" stats. Ignores "Q>=x" part.
        '''
        summary_dict = yaml.safe_load(f)

        # get q value threshold used for reads
        q_threshold = None
        for k in summary_dict.keys():
            if k.startswith('Q>='):
                q_threshold = k

        data_dict = {}
        data_dict['all'] = summary_dict['All reads']        # all reads
        data_dict['q_filt'] = summary_dict[q_threshold]     # quality filtered reads

        for q_key in ['all', 'q_filt']:
            for key_1 in ['reads', 'gigabases']:
                for key_2 in data_dict[q_key][key_1]:
                    new_key = '{} {}'.format(key_1, key_2)
                    data_dict[q_key][new_key] = data_dict[q_key][key_1][key_2]
                data_dict[q_key].pop(key_1)    # removes key after flattening

        return data_dict['all'], data_dict['q_filt'], q_threshold


    def headers_to_use(self):
        '''
        Defines features of columns to be used in multiqc table
        '''
        headers = OrderedDict()

        headers['total.reads'] = {
            'title': 'Total reads',
            'description': 'Total number of reads',
            'format': '{:,.0f}',
            'scale': 'Greys',
        }
        headers['total.gigabases'] = {
            'title': 'Total bases (GB)',
            'description': 'Total bases',
            'format': '{:,.2f}',
            'scale': 'Blues',
            # 'suffix': ' GB'
        }
        headers['N50.length'] = {
            'title': 'Reads N50',
            'description': 'Minimum read length needed to cover 50% of all reads',
            'format': '{:,.0f}',
            'scale': 'Purples',
        }
        headers['mean.q'] = {
            'title': 'Mean q',
            'description': 'Mean quality of reads',
            'min': 0,
            'max': 15,
            'format': '{:,.1f}',
            'hidden': True,
            'scale': 'Greens',
        }
        headers['median.q'] = {
            'title': 'Median q',
            'description': 'Median quality of reads',
            'min': 0,
            'max': 15,
            'format': '{:,.1f}',
            'scale': 'Greens',
        }
        headers['mean.length'] = {
            'title': 'Mean length',
            'description': 'Mean read length',
            'format': '{:,.0f}',
            'hidden': True,
            'scale': 'Blues',
        }
        headers['median.length'] = {
            'title': 'Median length',
            'description': 'Median read length',
            'format': '{:,.0f}',
            'scale': 'Blues',
            'hidden': True,
        }
        for s in ['>10kb', '>50kb', '>100kb']:
            hidden_status = False
            if s not in ['>10kb', '>50kb']:
                hidden_status = True

            headers['reads {}'.format(s)] = {
                'title': 'Total reads {}'.format(s),
                'description': 'Total number of reads {}'.format(s),
                'format': '{:,.0f}',
                'scale': 'Blues',
                'hidden': hidden_status
            }
        for s in ['>10kb', '>50kb', '>100kb']:
            headers['gigabases {}'.format(s)] = {
                'title': 'Bases {}'.format(s),
                'description': 'Total bases from reads {}'.format(s),
                'format': '{:,.3f}',
                'suffix': ' GB',
                'hidden': True,
                'scale': 'Blues',
            }

        return headers


    def table_qALL(self):
        """ Table showing stats for all reads """

        pconfig = { 'namespace': 'MinIONQC',
                        'id': 'minionqc-stats-qAll-table',
                        'table_title': 'minionqc-stats-qAll-table'}

        self.add_section (
            name = 'Stats: All reads',
            anchor = 'minionqc-stats-qAll',
            description = 'MinIONQC statistics for all reads',
            plot = table.plot(self.minionqc_data, self.headers_to_use(), pconfig=pconfig)
        )

        return None


    def table_qfiltered(self, q_threshold_list):
        """ Table showing stats for q-filtered reads """

        description = 'MinIONQC statistics for quality filtered reads. ' + \
                        'Quailty threshold used: {}.'.format(', '.join(list(q_threshold_list)))
        if len(q_threshold_list) > 1:
            description += '<br><b><i>Warning</b></i>: More than one quality thresholds were present.'
            log.warning('More than one quality thresholds were present. Thresholds: {}.'.format(', '.join(list(q_threshold_list))))

        pconfig = { 'namespace': 'MinIONQC',
                        'id': 'minionqc-stats-qFilt-table',
                        'table_title': 'eweffcref'}

        # 'rid' needs to be added to avoid linting error "HTML ID was a duplicate"
        headers = self.headers_to_use()
        for k in headers:
            headers[k]['rid'] = "rid_{}".format(headers[k]['title']) 

        self.add_section (
            name = 'Stats: Quality filtered reads',
            anchor = 'minionqc-stats-qFilt',
            description = description,
            plot = table.plot(self.qfilt_data, headers, pconfig=pconfig)
        )

        return None


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
