#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function

import json
import logging
import os
from collections import OrderedDict

import numpy

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.salmon.gcmodel import GCModel
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Salmon', anchor='salmon',
        href='http://combine-lab.github.io/salmon/',
        info="is a tool for quantifying the expression of transcripts using RNA-seq data.")

        # Parse meta information. JSON win!
        self.salmon_meta = dict()
        for f in self.find_log_files('salmon/meta'):
            # Get the s_name from the parent directory
            s_name = os.path.basename(os.path.dirname(f['root']))
            s_name = self.clean_s_name(s_name, f['root'])
            self.salmon_meta[s_name] = json.loads(f['f'])
        # Parse Fragment Length Distribution logs
        self.salmon_fld = dict()
        self.salmon_gc = []
        for f in self.find_log_files('salmon/fld'):
            # Get the s_name from the parent directory
            if os.path.basename(f['root']) == 'libParams':
                s_name = os.path.basename(os.path.dirname(f['root']))
                s_name = self.clean_s_name(s_name, f['root'])
                self.parse_gc_bias(f['root'])
                parsed = OrderedDict()
                for i, v in enumerate(f['f'].split()):
                    parsed[i] = float(v)
                if len(parsed) > 0:
                    if s_name in self.salmon_fld:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.salmon_fld[s_name] = parsed
        # Filter to strip out ignored sample names
        self.salmon_meta = self.ignore_samples(self.salmon_meta)
        self.salmon_fld = self.ignore_samples(self.salmon_fld)

        if len(self.salmon_meta) == 0 and len(self.salmon_fld) == 0:
            raise UserWarning
        if len(self.salmon_meta) > 0:
            log.info("Found {} meta reports".format(len(self.salmon_meta)))
            self.write_data_file(self.salmon_meta, 'multiqc_salmon')
        if len(self.salmon_fld) > 0:
            log.info("Found {} fragment length distributions".format(len(self.salmon_fld)))

        # Add alignment rate to the general stats table
        headers = OrderedDict()
        headers['percent_mapped'] = {
            'title': '% Aligned',
            'description': '% Mapped reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        headers['num_mapped'] = {
            'title': 'M Aligned',
            'description': 'Mapped reads (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: float(x) / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.salmon_meta, headers)

        # Fragment length distribution plot
        pconfig = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Fragment Length Distribution',
            'ylab': 'Fraction',
            'xlab': 'Fragment Length (bp)',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_fld, pconfig))
        self.plot_gc_bias()

    def plot_gc_bias(self):
        pconfig = {
            'smooth_points': 500,
            'id': 'salmon_gc_plot',
            'title': 'GC Bias',
            'ylab': 'Obs/Exp ratio',
            'xlab': 'bins',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        for sample_gc in self.salmon_gc:
            ratio = numpy.divide(sample_gc.obs_, sample_gc.exp_)
            plot_details = {
                0: OrderedDict(enumerate(ratio[0])),
                1: OrderedDict(enumerate(ratio[1])),
                2: OrderedDict(enumerate(ratio[2]))
            }
            self.add_section(plot=linegraph.plot(plot_details, pconfig))

    def parse_gc_bias(self, f_root):
        is_exp_gc_exists = os.path.exists(os.path.join(os.path.dirname(f_root), 'aux_info', 'exp_gc.gz'))
        is_obs_gc_exists = os.path.exists(os.path.join(os.path.dirname(f_root), 'aux_info', 'obs_gc.gz'))
        if is_exp_gc_exists and is_obs_gc_exists:
            gc = GCModel()
            gc.from_file(os.path.dirname(f_root))
            self.salmon_gc.append(gc)
