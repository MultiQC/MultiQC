#!/usr/bin/env python

""" MultiQC module to parse output from THetA2 """

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
        super(MultiqcModule, self).__init__(name='THetA2', anchor='theta2',
        href="http://compbio.cs.brown.edu/projects/theta/",
        info="<em>(Tumor Heterogeneity Analysis)</em> estimates tumour purity "\
             "and clonal / subclonal copy number.")

        # Find and load any THetA2 reports
        self.theta2_data = dict()
        for f in self.find_log_files('theta2', filehandles=True):
            parsed_data = self.parse_theta2_report(f['f'])
            if len(parsed_data) > 0:
                if f['s_name'] in self.theta2_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f)
                self.theta2_data[f['s_name']] = parsed_data

        # Filter to strip out ignored sample names
        self.theta2_data = self.ignore_samples(self.theta2_data)

        if len(self.theta2_data) == 0 :
            raise UserWarning

        log.info("Found {} reports".format(len(self.theta2_data)))

        # Write parsed report data to a file
        self.write_data_file(self.theta2_data, 'multiqc_theta2')

        # Alignment bar plot
        self.add_section (
            name = 'Tumour Subclone Purities',
            anchor = 'theta2-purities',
            description = 'Purities of tumour subclones. <em>NB:</em> Only first maximum likelihood solution for each sample shown.',
            plot = self.theta2_purities_chart()
        )


    def parse_theta2_report (self, fh):
        """ Parse the final THetA2 log file. """
        parsed_data = {}
        for l in fh:
            if l.startswith('#'):
                continue
            else:
                s = l.split("\t")
                purities = s[1].split(',')
                parsed_data['proportion_germline'] = float(purities[0]) * 100.0
                for i, v in enumerate(purities[1:]):
                    if i <= 5:
                        parsed_data['proportion_tumour_{}'.format(i+1)] = float(v) * 100.0
                    else:
                        parsed_data['proportion_tumour_gt5'] = (float(v) * 100.0) + parsed_data.get('proportion_tumour_gt5', 0)
                break
        return parsed_data


    def theta2_purities_chart (self):
        """ Make the plot showing alignment rates """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['proportion_germline'] =   { 'name': 'Germline' }
        keys['proportion_tumour_1'] =   { 'name': 'Tumour Subclone 1' }
        keys['proportion_tumour_2'] =   { 'name': 'Tumour Subclone 2' }
        keys['proportion_tumour_3'] =   { 'name': 'Tumour Subclone 3' }
        keys['proportion_tumour_4'] =   { 'name': 'Tumour Subclone 4' }
        keys['proportion_tumour_5'] =   { 'name': 'Tumour Subclone 5' }
        keys['proportion_tumour_gt5'] = { 'name': 'Tumour Subclones > 5' }

        # Config for the plot
        pconfig = {
            'id': 'theta2_purity_plot',
            'title': 'THetA2: Tumour Subclone Purities',
            'cpswitch': False,
            'ymin': 0,
            'ymax': 100,
            'ylab': '% Purity',
            'tt_suffix': '%'
        }

        return bargraph.plot(self.theta2_data, keys, pconfig)
