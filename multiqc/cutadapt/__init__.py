#!/usr/bin/env python

""" MultiQC module to parse output from Cutadapt """

from __future__ import print_function
import io
import logging
import os
import re

import multiqc
from multiqc import config

# Initialise the logger
log = logging.getLogger('MultiQC : {0:<14}'.format('Cutadapt'))

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(log)

        # Static variables
        self.name = "Cutadapt"
        self.anchor = "cutadapt"
        self.intro = '<p><a href="https://code.google.com/p/cutadapt/" target="_blank">Cutadapt</a> \
            is a tool to find and remove adapter sequences, primers, poly-A tails and other types \
            of unwanted sequence from your high-throughput sequencing reads.</p>'

        # Find and load any Cutadapt reports
        self.cutadapt_data = dict()
        self.cutadapt_length_counts = dict()
        self.cutadapt_length_exp = dict()
        self.cutadapt_length_obsexp = dict()
        
        for f in self.find_log_files(contents_match='This is cutadapt', filehandles=True):
            self.parse_cutadapt_logs(f)        

        if len(self.cutadapt_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.cutadapt_data)))

        # Write parsed report data to a file
        # Only the summary stats - skip the length data (t_lengths)
        with io.open (os.path.join(config.output_dir, 'report_data', 'multiqc_cutadapt.txt'), "w", encoding='utf-8') as f:
            print( self.dict_to_csv( self.cutadapt_data ), file=f)

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.cutadapt_general_stats_table()

        # Trimming Length Profiles
        # Only one section, so add to the intro
        self.intro += self.cutadapt_length_trimmed_plot()


    def parse_cutadapt_logs(self, f):
        """ Go through log file looking for cutadapt output """
        fh = f['f']
        regexes = {
            'bp_processed': "Total basepairs processed:\s*([\d,]+) bp",
            'bp_written': "Total written \(filtered\):\s*([\d,]+) bp",
            'quality_trimmed': "Quality-trimmed:\s*([\d,]+) bp",
            'r_processed': "Total reads processed:\s*([\d,]+)",
            'r_with_adapters': "Reads with adapters:\s*([\d,]+)"
        }
        s_name = None
        for l in fh:
            # New log starting
            if l.startswith('This is cutadapt'):
                s_name = None
            
            # Get sample name from end of command line params
            if l.startswith('Command line parameters'):
                s_name = l.split()[-1]
                s_name = self.clean_s_name(s_name, f['root'])
                if s_name in self.cutadapt_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.cutadapt_data[s_name] = dict()
                self.cutadapt_length_counts[s_name] = dict()
                self.cutadapt_length_exp[s_name] = dict()
                self.cutadapt_length_obsexp[s_name] = dict()
            
            if s_name is not None:
                # Search regexes for overview stats
                for k, r in regexes.items():
                    match = re.search(r, l)
                    if match:
                        self.cutadapt_data[s_name][k] = int(match.group(1).replace(',', ''))

                if 'length' in l and 'count' in l and 'expect' in l:
                    # Nested loop to read this section while the regex matches
                    for l in fh:
                        r_seqs = re.search("^(\d+)\s+(\d+)\s+([\d\.]+)", l)
                        if r_seqs:
                            a_len = int(r_seqs.group(1))
                            self.cutadapt_length_counts[s_name][a_len] = int(r_seqs.group(2))
                            self.cutadapt_length_exp[s_name][a_len] = float(r_seqs.group(3))
                            if float(r_seqs.group(3)) > 0:
                                self.cutadapt_length_obsexp[s_name][a_len] = float(r_seqs.group(2)) / float(r_seqs.group(3))
                            else:
                                # Cheating, I know. Infinity is difficult to plot.
                                self.cutadapt_length_obsexp[s_name][a_len] = float(r_seqs.group(2))
                        else:
                            break
        
        # Calculate a few extra numbers of our own
        for s_name in self.cutadapt_data.keys():
            if 'bp_processed' in self.cutadapt_data[s_name] and 'bp_written' in self.cutadapt_data[s_name]:
                self.cutadapt_data[s_name]['percent_trimmed'] = (float(self.cutadapt_data[s_name]['bp_processed'] - self.cutadapt_data[s_name]['bp_written']) / self.cutadapt_data[s_name]['bp_processed']) * 100



    def cutadapt_general_stats_table(self):
        """ Take the parsed stats from the Cutadapt report and add it to the
        basic stats table at the top of the report """

        config.general_stats['headers']['bp_trimmed'] = '<th class="chroma-col" data-chroma-scale="OrRd" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="% Total Base Pairs trimmed by Cutadapt">Trimmed</span></th>'
        for samp, vals in self.cutadapt_data.items():
            config.general_stats['rows'][samp]['bp_trimmed'] = '<td class="text-right">{:.1f}%</td>'.format(vals['percent_trimmed'])
    

    def cutadapt_length_trimmed_plot (self):
        """ Generate the trimming length plot """
        html = '<p>This plot shows the number of reads with certain lengths of adapter trimmed. \n\
        Obs/Exp shows the raw counts divided by the number expected due to sequencing errors. A defined peak \n\
        may be related to adapter length. See the \n\
        <a href="http://cutadapt.readthedocs.org/en/latest/guide.html#how-to-read-the-report" target="_blank">cutadapt documentation</a> \n\
        for more information on how these numbers are generated.</p>'
        
        pconfig = {
            'id': 'cutadapt_plot',
            'title': 'Lengths Trimmed',
            'ylab': 'Observed / Expected',
            'xlab': 'Length Trimmed (bp)',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [{'name': 'Obs/Exp', 'ylab': 'Observed / Expected'},
                            {'name': 'Counts', 'ylab': 'Count'}]
        }
        
        html += self.plot_xy_data([self.cutadapt_length_obsexp, self.cutadapt_length_counts], pconfig)
        
        return html
