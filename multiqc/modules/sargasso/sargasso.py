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

# xintodo main code goes here
class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='sargasso',
        anchor='sargasso', target='Subread featureCounts',
        href='http://statbio.github.io/Sargasso/',
        info="is a tool to separate mixed-species RNA-seq reads"\
             "according to their species of origin.")

         # Find and load any Sargasso reports
        self.sargasso_data = dict()
        self.sargasso_files = list()
        self.sargasso_keys = list() #header keys

        for f in self.find_log_files('sargasso'):
            self.parse_sargasso_logs(f)
            self.sargasso_files.append(f)
            # self.sargasso_data[f['s_name'].split(' | ')[-2]] = self.parse_sargasso_logs(f['f'])
            # self.sargasso_data[f['s_name']] = self.parse_sargasso_logs(f['f'])

        # log.info('Parsing log file...')
        # log.info(self.sargasso_data)

        #sample name has been changed, think again
        # log.info('Removing innored samples...')
        self.sargasso_data = self.ignore_samples(self.sargasso_data)

        if len(self.sargasso_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.sargasso_files)))

        # Write parsed report data to a file
        self.write_data_file(self.sargasso_data, 'multiqc_sargasso')

        # Basic Stats Table
        self.sargasso_stats_table()

        # Assignment bar plot
        self.add_section( plot = self.sargasso_chart() )

        # log.info('done')


    def _chunks(l, n):
        n = max(1, n)
        return (l[i:i+n] for i in xrange(0, len(l), n))

    def parse_sargasso_logs(self, f):
        """ Parse the sargasso log file. """
        # file_names = list()
        # log.info('Parsing ' + f['s_name'])
        parsed_data = dict()

        species_name = list()
        items = list()
        header = list()
        for l in f['f'].splitlines():
            thisrow = list()
            s = l.split(",")
            if len(s) < 7:
                continue
            if s[0] == 'Sample':
                header = s
                for i in header[1:]:
                    #find out what species included
                    sname = i.split('-')[-1]
                    if sname not in species_name:
                        species_name.append(sname)
                    #find out what is being counted
                    kname = ("-".join(i.split('-')[-3:-1]))
                    if kname not in items:
                        items.append(kname)
            else:
                #start sample lines.
                sample_name = s.pop(0)
                chunk_by_species = [s[i:i + len(items)] for i in xrange(0, len(s), len(items))];
                for idx,v in enumerate(chunk_by_species):
                    #adding species name to the same name for easy interpretation
                    new_sample_name = '_'.join([sample_name,species_name[idx]])

                    if new_sample_name in self.sargasso_data.keys():
                        log.debug("Duplicate sample name found! Overwriting: {}".format(new_sample_name))

                    try:
                        self.sargasso_data[new_sample_name] = dict(zip(items,map(int, v)))
                    except ValueError:
                        pass


        # Check that this actually is a Sargasso file
        if 'Assigned-Reads' not in items:
            return None

        self.sargasso_keys = items

        for idx, f_name in enumerate(self.sargasso_data.keys()):

            # Clean up sample name
            s_name = self.clean_s_name(f_name, f['root'])

            # Reorganised parsed data for this sample
            # Collect total read count number
            self.sargasso_data[f_name]['Total'] = 0;
            for key, value in self.sargasso_data[f_name].iteritems():   # iter on both keys and values
                if key.endswith("Reads"):
                    self.sargasso_data[f_name]['Total'] += value

            # Calculate the percent aligned if we can
            try:
                self.sargasso_data[f_name]['sargasso_percent_assigned'] = (float(self.sargasso_data[f_name]['Assigned-Reads'])/float(self.sargasso_data[f_name]['Total'])) * 100.0
            except (KeyError, ZeroDivisionError):
                pass


    def sargasso_stats_table(self):
        """ Take the parsed stats from the featureCounts report and add them to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['sargasso_percent_assigned'] = {
            'title': 'sargasso % Assigned',
            'description': 'sargasso % Assigned reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn'
        }
        headers['Assigned-Reads'] = {
            'title': 'sargasso {} Assigned'.format(config.read_count_prefix),
            'description': 'sargasso Assigned reads ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: float(x) * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.sargasso_data, headers)


    def sargasso_chart (self):
        """ Make the featureCounts assignment rates plot """

        # Config for the plot
        config = {
            'id': 'sargasso_assignment_plot',
            'title': 'Sargasso:',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        #We only want to plot the READs at the moment
        return bargraph.plot(self.sargasso_data, [name for name in self.sargasso_keys if 'Reads' in name], config)
