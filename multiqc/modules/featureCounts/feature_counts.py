#!/usr/bin/env python

""" MultiQC module to parse output from featureCounts """

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
        super(MultiqcModule, self).__init__(name='featureCounts',
        anchor='featurecounts', target='Subread featureCounts',
        href='http://bioinf.wehi.edu.au/featureCounts/',
        info="is a highly efficient general-purpose read summarization program"\
        " that counts mapped reads for genomic features such as genes, exons,"\
        " promoter, gene bodies, genomic bins and chromosomal locations.")

        # Find and load any featureCounts reports
        self.featurecounts_data = dict()
        self.featurecounts_keys = list()
        for f in self.find_log_files('featurecounts'):
            self.parse_featurecounts_report(f)

        # Filter to strip out ignored sample names
        self.featurecounts_data = self.ignore_samples(self.featurecounts_data)

        if len(self.featurecounts_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.featurecounts_data)))

        # Write parsed report data to a file
        self.write_data_file(self.featurecounts_data, 'multiqc_featureCounts')

        # Basic Stats Table
        self.featurecounts_stats_table()

        # Assignment bar plot
        self.add_section( plot = self.featureCounts_chart() )


    def parse_featurecounts_report (self, f):
        """ Parse the featureCounts log file. """

        file_names = list()
        parsed_data = dict()
        for l in f['f'].splitlines():
            thisrow = list()
            s = l.split("\t")
            if len(s) < 2:
                continue
            if s[0] == 'Status':
                for f_name in s[1:]:
                    file_names.append(f_name)
            else:
                k = s[0]
                if k not in self.featurecounts_keys:
                    self.featurecounts_keys.append(k)
                for val in s[1:]:
                    try:
                        thisrow.append(int(val))
                    except ValueError:
                        pass
            if len(thisrow) > 0:
                parsed_data[k] = thisrow
        # Check that this actually is a featureCounts file, as format and parsing is quite general
        if 'Assigned' not in parsed_data.keys():
            return None
        for idx, f_name in enumerate(file_names):

            # Clean up sample name
            s_name = self.clean_s_name(f_name, f['root'])

            # Reorganised parsed data for this sample
            # Collect total count number
            data = dict()
            data['Total'] = 0
            for k in parsed_data:
                data[k] = parsed_data[k][idx]
                data['Total'] += parsed_data[k][idx]

            # Calculate the percent aligned if we can
            try:
                data['percent_assigned'] = (float(data['Assigned'])/float(data['Total'])) * 100.0
            except (KeyError, ZeroDivisionError):
                pass

            # Add to the main dictionary
            if len(data) > 1:
                if s_name in self.featurecounts_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.add_data_source(f, s_name)
                self.featurecounts_data[s_name] = data


    def featurecounts_stats_table(self):
        """ Take the parsed stats from the featureCounts report and add them to the
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
        headers['Assigned'] = {
            'title': '{} Assigned'.format(config.read_count_prefix),
            'description': 'Assigned reads ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: float(x) * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.featurecounts_data, headers)


    def featureCounts_chart (self):
        """ Make the featureCounts assignment rates plot """

        # Config for the plot
        config = {
            'id': 'featureCounts_assignment_plot',
            'title': 'featureCounts Assignments',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(self.featurecounts_data, self.featurecounts_keys, config)
