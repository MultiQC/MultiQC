#!/usr/bin/env python

""" MultiQC module to parse output from HiCUP """

from __future__ import print_function
import logging

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='HiCUP', anchor='hicup',
        href='http://www.bioinformatics.babraham.ac.uk/projects/hicup/', 
        info="(Hi-C User Pipeline) is tool for mapping and performing "\
         "quality control on Hi-C data.")
        
        # Find and load any HiCUP summary reports
        self.hicup_data = dict()
        for f in self.find_log_files(config.sp['hicup']):
            self.parse_hicup_logs(f)

        if len(self.hicup_data) == 0:
            log.debug("Could not find any HiCUP data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.hicup_data)))
        
        import json
        print(json.dumps(self.hicup_data, indent=4))
        
        # Write parsed data to a file
        self.write_data_file(self.hicup_data, 'multiqc_hicup')

        # TODO: Make some sections and some plots


    def parse_hicup_logs(self, f):
        """ Parse a HiCUP summary report """
        if not f['fn'].endswith('.txt'):
            return None
        header = []
        lines = f['f'].splitlines()
        for l in lines:
            s = l.split("\t")
            if len(header) == 0:
                if s[0] != 'File':
                    return None
                header = s[1:]
            else:
                s_name = s[0].lstrip('HiCUP_output/')
                parsed_data = {}
                for idx, num in enumerate(s[1:]):
                    parsed_data[header[idx]] = num
                if s_name in self.hicup_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.add_data_source(f, s_name)
                self.hicup_data[s_name] = parsed_data


    def hicup_length_trimmed_plot (self):
        """ Generate the hicup plot """    
        pconfig = {
            'id': 'hicup_plot',
            'title': 'hicup complexity curve',
            'ylab': 'Unique Molecules',
            'xlab': 'Total Molecules (including duplicates)',
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
                 + self.plot_xy_data(self.hicup_data, pconfig)
