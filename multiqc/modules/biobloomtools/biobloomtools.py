#!/usr/bin/env python

""" MultiQC module to parse output from BioBloom Tools """

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
        super(MultiqcModule, self).__init__(name='BioBloom Tools', anchor='biobloomtools',
        href="https://github.com/bcgsc/biobloom/",
        info="creates filters for a genome reference and categorises "\
             "sequences. This is faster than alignment and can be used for pre-processing "\
             "and QC applications such as contamination detection.")

        # Find and load any BioBloom Tools reports
        self.bbt_data = dict()
        self.num_orgs = 0
        for f in self.find_log_files('biobloomtools', filehandles=True):
            parsed_data = self.parse_bbt(f['f'])
            if len(parsed_data) > 0:
                if f['s_name'] in self.bbt_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f)
                self.bbt_data[f['s_name']] = parsed_data

        # Filter to strip out ignored sample names
        self.bbt_data = self.ignore_samples(self.bbt_data)

        if len(self.bbt_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.bbt_data)))

        # Section 1 - Alignment Profiles
        self.add_section( plot = self.bbt_simple_plot() )

        # Write the total counts and percentages to files
        data_export = {}
        for s_name in self.bbt_data:
            data_export[s_name] = dict()
            for org in self.bbt_data[s_name]:
                data_export[s_name][org] = int(self.bbt_data[s_name][org]['hits'])
        self.write_data_file(data_export, 'multiqc_biobloomtools')


    def parse_bbt(self, fh):
        """ Parse the BioBloom Tools output into a 3D dict """
        parsed_data = OrderedDict()
        headers = None
        for l in fh:
            s = l.split("\t")
            if headers is None:
                headers = s
            else:
                parsed_data[s[0]] = dict()
                for i, h in enumerate(headers[1:]):
                    parsed_data[s[0]][h] = float(s[i+1])

        return parsed_data

    def bbt_simple_plot(self):
        """ Makes a simple bar plot with summed alignment counts for
        each species, stacked. """

        # First, sum the different types of alignment counts
        data = OrderedDict()
        cats = OrderedDict()
        for s_name in self.bbt_data:
            data[s_name] = OrderedDict()
            for org in self.bbt_data[s_name]:
                data[s_name][org] = self.bbt_data[s_name][org]['hits'] - self.bbt_data[s_name][org]['shared']
                if org not in cats and org != 'multiMatch' and org != 'noMatch':
                    if org.lower().endswith('.fa'):
                        cname = org[:-3]
                    elif org.lower().endswith('.fasta'):
                        cname = org[:-6]
                    else:
                        cname = org
                    cats[org] = { 'name': cname }

        pconfig = {
            'id': 'biobloom_tools',
            'title': 'BioBloom Tools: Alignment counts per species',
            'ylab': 'Number of hits',
            'hide_zero_cats': False
        }
        cats['multiMatch'] = { 'name': 'Multiple Genomes', 'color': '#820000' }
        cats['noMatch'] = { 'name': 'No Match', 'color': '#cccccc' }

        return bargraph.plot(data, cats, pconfig)
