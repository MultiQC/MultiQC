#!/usr/bin/env python

""" MultiQC module to parse output from Kallisto """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Kallisto', anchor='kallisto',
        href="http://pachterlab.github.io/kallisto/",
        info="is a program for quantifying abundances of transcripts from RNA-Seq data.")

        # Find and load any Kallisto reports
        self.kallisto_data = dict()
        for f in self.find_log_files(config.sp['kallisto'], filehandles=True):
            self.parse_kallisto_log(f)

        if len(self.kallisto_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.kallisto_data)))

        # Write parsed report data to a file
        self.write_data_file(self.kallisto_data, 'multiqc_kallisto')

        # Basic Stats Table
        self.kallisto_general_stats_table()

        # Alignment Rate Plot
        self.intro += self.kallisto_alignment_plot()


    def parse_kallisto_log(self, f):
        s_name = total_reads = paligned_reads = fraglength = None
        for l in f['f']:
            
            # Get input filename
            match = re.search(r'\[quant\] will process pair 1: (\S+)', l)
            if match:
                s_name = self.clean_s_name(match.group(1), f['root'])
            
            if s_name is not None:
                aligned = re.search(r'\[quant\] processed ([\d,]+) reads, ([\d,]+) reads pseudoaligned', l)
                if aligned:
                    total_reads = float(aligned.group(1).replace(',',''))
                    paligned_reads = float(aligned.group(2).replace(',',''))
                flength = re.search(r'\[quant\] estimated average fragment length: ([\d\.]+)', l)
                if flength:
                    fraglength = float(flength.group(1).replace(',',''))
                
                if total_reads is not None and paligned_reads is not None and fraglength is not None:
                    if s_name in self.kallisto_data:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.kallisto_data[s_name] = {
                        'total_reads': total_reads,
                        'pseudoaligned_reads': paligned_reads,
                        'not_pseudoaligned_reads': total_reads - paligned_reads,
                        'percent_aligned': (paligned_reads / total_reads) * 100,
                        'fragment_length': fraglength,
                    }
                    s_name = total_reads = paligned_reads = fraglength = None

    def kallisto_general_stats_table(self):
        """ Take the parsed stats from the Kallisto report and add it to the
        basic stats table at the top of the report """
        
        headers = OrderedDict()
        headers['fragment_length'] = {
            'title': 'Frag Length',
            'description': 'Estimated average fragment length',
            'min': 0,
            'suffix': 'bp',
            'scale': 'RdYlGn',
            'format': '{:.1f}',
        }
        headers['percent_aligned'] = {
            'title': '% Aligned',
            'description': '% processed reads that were pseudoaligned',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
            'format': '{:.1f}%'
        }
        headers['pseudoaligned_reads'] = {
            'title': 'M Aligned',
            'description': 'pseudoaligned reads (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.kallisto_data, headers)

    def kallisto_alignment_plot (self):
        """ Make the HighCharts HTML to plot the alignment rates """
        
        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['pseudoaligned_reads'] = { 'color': '#437bb1', 'name': 'Pseudoaligned' }
        keys['not_pseudoaligned_reads'] =   { 'color': '#b1084c', 'name': 'Not aligned' }
        
        # Config for the plot
        config = {
            'title': 'Kallisto Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return plots.bargraph.plot(self.kallisto_data, keys, config)
