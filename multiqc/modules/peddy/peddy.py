#!/usr/bin/env python

""" MultiQC module to parse output from Peddy """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Peddy module class, parses stderr logs.
    Also understands logs saved by Trim Galore!
    (which contain peddy logs)
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Peddy', anchor='peddy',
        href='https://github.com/brentp/peddy',
        info="calculates genotype :: pedigree correspondence checks, ancestry checks and sex checks using VCF files.")

        # Find and load any Peddy reports
        self.peddy_data = dict()
        self.peddy_length_counts = dict()
        self.peddy_length_exp = dict()
        self.peddy_length_obsexp = dict()

        for f in self.find_log_files(config.sp['peddy']['summary_table']):
            parsed_data = self.parse_peddy_summary(f)
            if parsed_data is not None:
                self.peddy_data.update(parsed_data)

        if len(self.peddy_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.peddy_data)))

        # Write parsed report data to a file
        self.write_data_file(self.peddy_data, 'multiqc_peddy')

        # Basic Stats Table
        self.peddy_general_stats_table()

        self.sections = list()

        # PCA plot
        self.sections.append({
            'name': 'PCA Projection onto 1000 Genomes',
            'anchor': 'peddy-pca-plot',
            'content': self.peddy_pca_plot()
        })


    def parse_peddy_summary(self, f):
        """ Go through log file looking for peddy output """
        parsed_data = dict()
        headers = None
        for l in f['f'].splitlines():
            s = l.split("\t")
            if headers is None:
                s[0] = s[0].lstrip('#')
                headers = s
            else:
                parsed_data[s[1]] = dict()
                for i, v in enumerate(s):
                    if i != 1:
                        parsed_data[s[1]][headers[i]] = v
        if len(parsed_data) == 0:
            return None
        return parsed_data


    def peddy_general_stats_table(self):
        """ Take the parsed stats from the Peddy report and add it to the
        basic stats table at the top of the report """

        family_ids = [ x['family_id'] for x in self.peddy_data.values() ]

        headers = OrderedDict()
        headers['family_id'] = {
            'title': 'Family ID',
            'hidden': True if all([v == family_ids[0] for v in family_ids]) else False
        }
        headers['ancestry-prediction'] = {
            'title': 'Ancestry',
            'description': 'Ancestry Prediction',
        }
        headers['sex_het_ratio'] = {
            'title': 'Sex / Het Ratio',
        }
        headers['sex_error'] = {
            'title': 'Sex Error',
            'description': 'Error in sample sex prediction',
        }
        self.general_stats_addcols(self.peddy_data, headers)

    def peddy_pca_plot(self):

        data = dict()
        for s_name, d in self.peddy_data.items():
            data[s_name] = {
                'x': float(d['PC1']),
                'y': float(d['PC2']),
            }

        return plots.scatter.plot(data)

