#!/usr/bin/env python

""" MultiQC module to parse output from Preseq """

from __future__ import print_function
import logging

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Preseq', anchor='preseq',
        href='http://smithlabresearch.org/software/preseq/', 
        info="estimates the complexity of a library, showing how many additional "\
         "unique reads are sequenced for increasing total read count.")
        
        # Find and load any Preseq reports
        self.preseq_data = dict()
        self.total_max = 0
        for f in self.find_log_files(config.sp['preseq']):
            parsed_data = self.parse_preseq_logs(f)
            if parsed_data is not None:
                if f['s_name'] in self.preseq_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f)
                self.preseq_data[f['s_name']] = parsed_data

        if len(self.preseq_data) == 0:
            log.debug("Could not find any preseq data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.preseq_data)))

        # Preseq plot
        # Only one section, so add to the intro
        self.intro += self.preseq_length_trimmed_plot()


    def parse_preseq_logs(self, f):
        """ Go through log file looking for preseq output """
        
        lines = f['f'].splitlines()
        header = lines.pop(0)
        if header.startswith('TOTAL_READS	EXPECTED_DISTINCT'):
            self.axis_label = 'Molecules'
        elif header.startswith('TOTAL_BASES	EXPECTED_DISTINCT'):
            self.axis_label = 'Bases'
        else:
            log.debug("First line of preseq file {} did not look right".format(f['fn']))
            return None
        
        data = {}
        for l in lines:
            s = l.split()
            # Sometimes the Expected_distinct count drops to 0, not helpful
            if float(s[1]) == 0 and float(s[0]) > 0:
                continue
            data[float(s[0])] = float(s[1])
            self.total_max = max(float(s[1]), self.total_max)
        return data


    def preseq_length_trimmed_plot (self):
        """ Generate the preseq plot """    
        pconfig = {
            'id': 'preseq_plot',
            'title': 'Preseq complexity curve',
            'ylab': 'Unique {}'.format(self.axis_label),
            'xlab': 'Total {} (including duplicates)'.format(self.axis_label),
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} total</b>: {point.y:,.0f} unique',
            'extra_series': [{
                'name': 'x = y',
                'data': [[0, 0], [self.total_max, self.total_max]],
                'dashStyle': 'Dash',
                'lineWidth': 1,
                'color': '#000000',
                'marker': { 'enabled': False },
                'enableMouseTracking': False,
                'showInLegend': False,
            }]
        }
        return "<p>A shallow curve indicates complexity saturation. The dashed line \
                shows a perfectly complex library where total reads = unique reads.</o>" \
                 + plots.linegraph.plot(self.preseq_data, pconfig)
