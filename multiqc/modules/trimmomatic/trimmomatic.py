#!/usr/bin/env python

""" MultiQC module to parse output from Trimmomatic """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Trimmomatic', anchor='trimmomatic',
        href='http://www.usadellab.org/cms/?page=trimmomatic',
        info="is a flexible read trimming tool for Illumina NGS data.")

        # Parse logs
        self.trimmomatic = dict()
        for f in self.find_log_files('trimmomatic', filehandles=True):
            self.parse_trimmomatic(f)

        # Filter to strip out ignored sample names
        self.trimmomatic = self.ignore_samples(self.trimmomatic)

        if len(self.trimmomatic) == 0:
            log.debug("Could not find any Trimmomatic data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} logs".format(len(self.trimmomatic)))
        self.write_data_file(self.trimmomatic, 'multiqc_trimmomatic')

        # Add drop rate to the general stats table
        headers = OrderedDict()
        headers['dropped_pct'] = {
            'title': '% Dropped',
            'description': '% Dropped reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'OrRd'
        }
        self.general_stats_addcols(self.trimmomatic, headers)

        # Make barplot
        self.trimmomatic_barplot()

    def parse_trimmomatic(self, f):
        s_name = None
        for l in f['f']:
            # Get the sample name
            if 'Trimmomatic' in l and 'Started with arguments:' in l:
                # Match everything up until the first .fastq or .fq
                match = re.search('Trimmomatic[SP]E: Started with arguments:.+?(?=\.fastq|\.fq)', l)
                if match:
                    # backtrack from the end to the first space
                    s_name = match.group().split()[-1]
                    s_name = self.clean_s_name(s_name, f['root'])
                    if s_name in self.trimmomatic:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                else:
                    # Try looking on the next line instead, sometimes have a line break (see issue #212)
                    l = next(f['f'])
                    match = re.search('.+?(?=\.fastq|\.fq)', l)
                    if match:
                        # backtrack from the end to the first space
                        s_name = match.group().split()[-1]
                        s_name = self.clean_s_name(s_name, f['root'])
                        if s_name in self.trimmomatic:
                            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))

            # Get single end stats
            if 'Input Reads' in l and s_name is not None:
                match = re.search('Input Reads: (\d+) Surviving: (\d+) \(([\d\.,]+)%\) Dropped: (\d+) \(([\d\.,]+)%\)', l)
                if match:
                    self.trimmomatic[s_name] = dict()
                    self.trimmomatic[s_name]['input_reads'] = float( match.group(1) )
                    self.trimmomatic[s_name]['surviving'] = float( match.group(2) )
                    self.trimmomatic[s_name]['surviving_pct'] = float( match.group(3).replace(',','.') )
                    self.trimmomatic[s_name]['dropped'] = float( match.group(4) )
                    self.trimmomatic[s_name]['dropped_pct'] = float( match.group(5).replace(',','.') )
                    s_name = None

            # Get paired end stats
            if 'Input Read Pairs' in l and s_name is not None:
                match = re.search('Input Read Pairs: (\d+) Both Surviving: (\d+) \(([\d\.,]+)%\) Forward Only Surviving: (\d+) \(([\d\.,]+)%\) Reverse Only Surviving: (\d+) \(([\d\.,]+)%\) Dropped: (\d+) \(([\d\.,]+)%\)', l)
                if match:
                    self.trimmomatic[s_name] = dict()
                    self.trimmomatic[s_name]['input_read_pairs'] = float( match.group(1) )
                    self.trimmomatic[s_name]['surviving'] = float( match.group(2) )
                    self.trimmomatic[s_name]['surviving_pct'] = float( match.group(3).replace(',','.') )
                    self.trimmomatic[s_name]['forward_only_surviving'] = float( match.group(4) )
                    self.trimmomatic[s_name]['forward_only_surviving_pct'] = float( match.group(5).replace(',','.') )
                    self.trimmomatic[s_name]['reverse_only_surviving'] = float( match.group(6) )
                    self.trimmomatic[s_name]['reverse_only_surviving_pct'] = float( match.group(7).replace(',','.') )
                    self.trimmomatic[s_name]['dropped'] = float( match.group(8) )
                    self.trimmomatic[s_name]['dropped_pct'] = float( match.group(9).replace(',','.') )
                    s_name = None


    def trimmomatic_barplot (self):
        """ Make the HighCharts HTML to plot the trimmomatic rates """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['surviving'] =              { 'color': '#437bb1', 'name': 'Surviving Reads' }
        keys['both_surviving'] =         { 'color': '#f7a35c', 'name': 'Both Surviving' }
        keys['forward_only_surviving'] = { 'color': '#e63491', 'name': 'Forward Only Surviving' }
        keys['reverse_only_surviving'] = { 'color': '#b1084c', 'name': 'Reverse Only Surviving' }
        keys['dropped'] =                { 'color': '#7f0000', 'name': 'Dropped' }

        # Config for the plot
        pconfig = {
            'id': 'trimmomatic_plot',
            'title': 'Trimmomatic',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        self.add_section( plot = bargraph.plot(self.trimmomatic, keys, pconfig) )
