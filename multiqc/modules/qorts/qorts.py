#!/usr/bin/env python

""" MultiQC module to parse output from QoRTs """

from __future__ import print_function
from collections import OrderedDict
import re
import logging

from multiqc import config
from multiqc.plots import bargraph, beeswarm
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='QoRTs', anchor='qorts',
        href='http://hartleys.github.io/QoRTs/',
        info="is toolkit for analysis, QC and data management of RNA-Seq datasets.")

        # Parse logs
        self.qorts_data = dict()
        for f in self.find_log_files('qorts', filehandles=True):
            self.parse_qorts(f)

        # Remove empty samples
        self.qorts_data = { s:v for s, v in self.qorts_data.items() if len(v) > 0 }

        # Filter to strip out ignored sample names
        self.qorts_data = self.ignore_samples(self.qorts_data)

        if len(self.qorts_data) == 0:
            log.debug("Could not find any QoRTs data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} logs".format(len(self.qorts_data)))
        self.write_data_file(self.qorts_data, 'multiqc_qorts')

        # Make plots
        self.qorts_splice_loci_barplot()
        self.qorts_splice_events_barplot()
        self.qorts_genebodycoverage_plot()

    def parse_qorts(self, f):
        s_names = None
        for l in f['f']:
            s = l.split("\t")
            if s_names is None:
                s_names = [ self.clean_s_name(s_name, f['root']) for s_name in s[1:] ]
                for s_name in s_names:
                    if s_name in self.qorts_data:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.qorts_data[s_name] = dict()
            else:
                for i, s_name in enumerate(s_names):
                    self.qorts_data[s_name][s[0]] = float(s[i+1])

    def qorts_splice_loci_barplot (self):
        """ Make the HighCharts HTML to plot the qorts splice loci """
        # Specify the order of the different possible categories
        keys = [
            'SpliceLoci_Known_ManyReads',
            'SpliceLoci_Known_FewReads',
            'SpliceLoci_Known_NoReads',
            'SpliceLoci_Novel_ManyReads',
            'SpliceLoci_Novel_FewReads',
        ]
        cats = OrderedDict()
        for k in keys:
            name = k.replace('SpliceLoci_', '').replace('_',': ')
            name = re.sub("([a-z])([A-Z])","\g<1> \g<2>",name)
            cats[k] = { 'name': name }

        # Config for the plot
        pconfig = {
            'id': 'qorts_splice_loci',
            'title': 'QoRTs: Splice Loci',
            'ylab': '# Splice Loci',
            'cpswitch_counts_label': 'Number of Splice Loci',
            'hide_zero_cats': False
        }

        self.add_section(
            name = "Splice Loci",
            plot = bargraph.plot(self.qorts_data, cats, pconfig)
        )

    def qorts_splice_events_barplot (self):
        """ Make the HighCharts HTML to plot the qorts splice events """
        # Specify the order of the different possible categories
        keys = [
            'SpliceEvents_KnownLociWithManyReads',
            'SpliceEvents_KnownLociWithFewReads',
            'SpliceEvents_NovelLociWithManyReads',
            'SpliceEvents_NovelLociWithFewReads',
        ]
        cats = OrderedDict()
        for k in keys:
            name = k.replace('SpliceEvents_', '')
            name = re.sub("([a-z])([A-Z])","\g<1> \g<2>",name)
            cats[k] = { 'name': name }

        # Config for the plot
        pconfig = {
            'id': 'qorts_splice_events',
            'title': 'QoRTs: Splice Events',
            'ylab': '# Splice Events',
            'cpswitch_counts_label': 'Number of Splice Events',
            'hide_zero_cats': False
        }

        self.add_section(
            name = "Splice Events",
            plot = bargraph.plot(self.qorts_data, cats, pconfig)
        )

    def qorts_genebodycoverage_plot (self):
        """ Make a beeswarm plot of the GeneBodyCoverage values """

        keys = [
            'GeneBodyCoverage_Overall_Mean',
            'GeneBodyCoverage_Overall_Median',
            'GeneBodyCoverage_LowExpress_Mean',
            'GeneBodyCoverage_LowExpress_Median',
            'GeneBodyCoverage_UMQuartile_Mean',
            'GeneBodyCoverage_UMQuartile_Median'
        ]
        cats = OrderedDict()
        for k in keys:
            name = k.replace('GeneBodyCoverage_', '')
            name = name.replace('_', ' ')
            name = re.sub("([a-z])([A-Z])","\g<1> \g<2>",name)
            cats[k] = {
                'title': name,
                'min': 0,
                'max': 1,
            }

        # Config for the plot
        pconfig = {
            'id': 'qorts_gene_body_coverage',
            'title': 'QoRTs: Gene Body Coverage'
        }

        self.add_section(
            name = 'Gene Body Coverage',
            plot = beeswarm.plot(self.qorts_data, cats, pconfig)
        )
