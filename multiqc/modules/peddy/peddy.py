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
                for s_name in parsed_data:
                    try:
                        self.peddy_data[s_name].update(parsed_data[s_name])
                    except KeyError:
                        self.peddy_data[s_name] = parsed_data[s_name]

        for f in self.find_log_files(config.sp['peddy']['sex_check']):
            parsed_data = self.parse_peddy_sexcheck(f)
            if parsed_data is not None:
                for s_name in parsed_data:
                    try:
                        self.peddy_data[s_name].update(parsed_data[s_name])
                    except KeyError:
                        self.peddy_data[s_name] = parsed_data[s_name]

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
        pca_plot = self.peddy_pca_plot()
        if pca_plot is not None:
            self.sections.append({
                'name': 'PCA Plot',
                'anchor': 'peddy-pca-plot',
                'content': pca_plot
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

    def parse_peddy_sexcheck(self, f):
        """ Go through the sex check output from peddy """
        parsed_data = dict()
        headers = None
        s_name_idx = None
        for l in f['f'].splitlines():
            s = l.split(",")
            if headers is None:
                headers = s
                s_name_idx = headers.index("sample_id")
            else:
                parsed_data[s[s_name_idx]] = dict()
                for i, v in enumerate(s):
                    if i != s_name_idx:
                        parsed_data[s[s_name_idx]][headers[i]] = v
        if len(parsed_data) == 0:
            return None
        return parsed_data

    def peddy_general_stats_table(self):
        """ Take the parsed stats from the Peddy report and add it to the
        basic stats table at the top of the report """

        family_ids = [ x.get('family_id') for x in self.peddy_data.values() ]

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
        headers['error'] = {
            'title': 'Sex Error',
            'description': 'Error in sample sex prediction',
        }
        self.general_stats_addcols(self.peddy_data, headers)

    def peddy_pca_plot(self):

        data = dict()
        for s_name, d in self.peddy_data.items():
            if 'PC1' in d and 'PC2' in d:
                data[s_name] = {
                    'x': float(d['PC1']),
                    'y': float(d['PC2']),
                }

        pconfig = {
            'id': 'peddy_pca_plot',
            'title': 'Peddy PCA Plot',
        }

        if len(data) > 0:
            return plots.scatter.plot(data, pconfig)

