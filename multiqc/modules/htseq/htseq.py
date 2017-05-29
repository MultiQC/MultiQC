#!/usr/bin/env python

""" MultiQC module to parse output from HTSeq Count """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='HTSeq Count',
        anchor='htseq', target='HTSeq Count',
        href='http://www-huber.embl.de/HTSeq/doc/count.html',
        info=" is part of the HTSeq Python package - it takes a file with aligned sequencing "\
             "reads, plus a list of genomic features and counts how many reads map to each feature.")

        # Find and load any HTSeq Count reports
        self.htseq_data = dict()
        self.htseq_keys = list()
        for f in self.find_log_files('htseq', filehandles=True):
            parsed_data = self.parse_htseq_report(f)
            if parsed_data is not None:
                self.htseq_data[f['s_name']] = parsed_data

        # Filter to strip out ignored sample names
        self.htseq_data = self.ignore_samples(self.htseq_data)

        if len(self.htseq_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.htseq_data)))

        # Write parsed report data to a file
        self.write_data_file(self.htseq_data, 'multiqc_htseq')

        # Basic Stats Table
        self.htseq_stats_table()

        # Assignment bar plot
        self.add_section( plot = self.htseq_counts_chart() )


    def parse_htseq_report (self, f):
        """ Parse the HTSeq Count log file. """
        keys = [ '__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique' ]
        parsed_data = dict()
        assigned_counts = 0
        for l in f['f']:
            s = l.split("\t")
            if s[0] in keys:
                parsed_data[s[0][2:]] = int(s[1])
            else:
                try:
                    assigned_counts += int(s[1])
                except (ValueError, IndexError):
                    pass
        if len(parsed_data) > 0:
            parsed_data['assigned'] = assigned_counts
            parsed_data['total_count'] = sum([v for v in parsed_data.values()])
            parsed_data['percent_assigned'] = (float(parsed_data['assigned']) / float(parsed_data['total_count'])) * 100.0
            return parsed_data
        return None


    def htseq_stats_table(self):
        """ Take the parsed stats from the HTSeq Count report and add them to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['percent_assigned'] = {
            'title': '% Assigned',
            'description': '% Assigned reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn'
        }
        headers['assigned'] = {
            'title': '{} Assigned'.format(config.read_count_prefix),
            'description': 'Assigned Reads ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: float(x) * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.htseq_data, headers)


    def htseq_counts_chart (self):
        """ Make the HTSeq Count assignment rates plot """
        cats = OrderedDict()
        cats['assigned'] =      { 'name': 'Assigned' }
        cats['ambiguous'] =     { 'name': 'Ambiguous' }
        cats['alignment_not_unique'] = { 'name': 'Alignment Not Unique' }
        cats['no_feature'] =    { 'name': 'No Feature' }
        cats['too_low_aQual'] = { 'name': 'Too Low aQual' }
        cats['not_aligned'] =   { 'name': 'Not Aligned' }
        config = {
            'id': 'htseq_assignment_plot',
            'title': 'HTSeq Count Assignments',
            'ylab': '# Reads',
            'hide_zero_cats': False,
            'cpswitch_counts_label': 'Number of Reads'
        }
        return bargraph.plot(self.htseq_data, cats, config)
