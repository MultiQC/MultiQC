#!/usr/bin/env python

""" MultiQC module to parse output from ngsderive """

from __future__ import print_function
import csv
import io
import logging
import re

from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    ngsderive module class, parses stderr logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='ngsderive', anchor='ngsderive',
        href='https://github.com/claymcleod/ngsderive',
        info="is a forensic analysis tool useful in backwards computing information from " + \
             "next-generation sequencing data. " + \
             "Notably, results are provided as a 'best guess' — the tool does " + \
             "not claim 100%% accuracy and results should be considered with that understanding. " + \
             "Please see the <a href='https://github.com/claymcleod/ngsderive/#ngsderive'>documentation</a> " + \
             "for more information. ")

        self.strandedness = {}
        self.instrument = {}
        self.readlen = {}

        # parse ngsderive summary file
        for f in self.find_log_files('ngsderive/strandedness'):
            self.parse(self.strandedness, f.get('f'), f.get('s_name'))

        for f in self.find_log_files('ngsderive/instrument'):
            self.parse(self.instrument, f.get('f'), f.get('s_name'))

        for f in self.find_log_files('ngsderive/readlen'):
            self.parse(self.readlen, f.get('f'), f.get('s_name'))

        self.strandedness = self.ignore_samples(self.strandedness)
        self.instrument = self.ignore_samples(self.instrument)
        self.readlen = self.ignore_samples(self.readlen)

        any_results_found = False
        if self.strandedness:
            self.add_strandedness_data()
            any_results_found = True

        if self.instrument:
            self.add_instrument_data()
            any_results_found = True

        if self.readlen:
            self.add_readlen_data()
            any_results_found = True

        if not any_results_found:
            raise UserWarning


    def parse(self, d, contents, sample_name):
        relevant_items = []

        f = io.StringIO(contents)
        dialect = csv.Sniffer().sniff(f.read(1024))
        f.seek(0)

        for row in csv.DictReader(io.StringIO(contents), dialect=dialect):
            if sample_name in row.get("File"):
                relevant_items.append(row)

        if len(relevant_items) < 1:
            raise RuntimeError("Could not find results for {sample_name} in ngsderive strandedness report!".format(sample_name=sample_name))
        elif len(relevant_items) > 1:
            raise RuntimeError("Too many results for {sample_name} in ngsderive strandedness report!".format(sample_name=sample_name))

        d[sample_name] = relevant_items.pop(0)


    def add_strandedness_data(self):
        data = {}
        for sample, strandedness in self.strandedness.items():
            data[sample] = {
                "predicted": strandedness.get("Predicted"),
                "forward": round(float(strandedness.get("ForwardPct")) * 100.0, 2),
                "reverse": round(float(strandedness.get("ReversePct")) * 100.0, 2),
            }

        headers = OrderedDict()
        headers['predicted'] = {
            'title': 'Strandedness (Predicted)',
            'description': 'Predicted strandedness from ngsderive'
        }
        headers['forward'] = {
            'title': 'Strandedness (Forward %)',
            'description': 'Percentage of reads which were evidence for forward-stranded.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'RdYlGn',
            'hidden': True
        }
        headers['reverse'] = {
            'title': 'Strandedness (Reverse %)',
            'description': 'Percentage of reads which were evidence for reverse-stranded.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'RdYlGn',
            'hidden': True
        }
        self.general_stats_addcols(data, headers)

        # Config for the plot
        pconfig = {
            'id': 'ngsderive_strandedness_plot',
            'title': 'ngsderive: strandedness',
            'ylab': '% Read Evidence',
            'ymin': 0,
            'ymax': 100,
            'tt_percentages': True,
            'ylab_format': '{value}%',
            'cpswitch': False
        }

        self.add_section (
            name = 'Strandedness',
            anchor = 'ngsderive-strandedness',
            description = 'Predicted strandedness provided by ngsderive. ' + \
                'For more information, please see <a href="https://github.com/claymcleod/ngsderive/#strandedness-inference">' + \
                'the relevant documentation and limitations</a>.',
            plot = bargraph.plot(data, headers, pconfig)
        )

    def add_instrument_data(self):
        data = {}
        for sample, instrument in self.instrument.items():
            data[sample] = {
                "instrument": '/'.join(sorted(instrument.get("Instrument").split(" or "))),
                "confidence": instrument.get("Confidence"),
                "basis": instrument.get("Basis"),
            }

        headers = OrderedDict()
        headers['instrument'] = {
            'title': 'Instrument (Predicted)',
            'description': 'Predicted instrument from ngsderive'
        }
        headers['confidence'] = {
            'title': 'Instrument (Confidence)',
            'description': 'Level of confidence (low, medium, high) that the predicted instrument is correct.',
        }
        headers['basis'] = {
            'title': 'Instrument (Basis)',
            'description': 'Basis upon which the prediction was made (see documentation for more details).',
        }
        self.general_stats_addcols(data, headers)


    def add_readlen_data(self):
        data = {}
        for sample, readlen in self.readlen.items():
            data[sample] = {
                "evidence": readlen.get("Evidence"),
                "majoritypctdetected": round(float(readlen.get("MajorityPctDetected")) * 100.0, 2),
                "consensusreadlength": int(readlen.get("ConsensusReadLength")),
            }

        headers = OrderedDict()
        headers['consensusreadlength'] = {
            'title': 'Read Length (Predicted)',
            'description': 'Predicted read length from ngsderive.',
            'scale': False,
            'format': '{:,.d}',
        }
        headers['majoritypctdetected'] = {
            'title': 'Read Length (% Supporting Reads)',
            'description': 'Percentage of reads which were measured at the predicted read length.',
        }
        headers['evidence'] = {
            'title': 'Read Length (Evidence)',
            'description': 'Evidence from readlen analysis (see documentation for more details). ' + \
                           'Typically, this is not helpful to the user in text form.',
            'hidden': True
        }
        self.general_stats_addcols(data, headers)

        linedata = {}
        for sample, d in data.items():
            linedata[sample] = {}
            for parts in d.get("evidence").split(";"):
                (k, v) = parts.split("=")
                linedata[sample][int(k)] = int(v)

        # Config for the plot
        pconfig = {
            'id': 'ngsderive_readlen_plot',
            'title': 'ngsderive: readlen',
            'ylab': 'Evidence for Respective Read Length'
        }

        self.add_section (
            name = 'Read length',
            anchor = 'ngsderive-readlen',
            description = 'Predicted read length provided by ngsderive. ' + \
                'For more information, please see <a href="https://github.com/claymcleod/ngsderive/#read-length-calculation">' + \
                'the relevant documentation and limitations</a>.',
            plot = linegraph.plot(linedata, pconfig=pconfig)
        )