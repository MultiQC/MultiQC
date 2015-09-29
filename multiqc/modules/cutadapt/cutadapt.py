#!/usr/bin/env python

""" MultiQC module to parse output from Cutadapt """

from __future__ import print_function
import io
import logging
import os
import re

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger()

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Cutadapt', anchor='cutadapt',
        href='https://code.google.com/p/cutadapt/', 
        info="is a tool to find and remove adapter sequences, primers, poly-A"\
         "tails and other types of unwanted sequence from your high-throughput"\
         " sequencing reads.")

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

        headers = {}
        headers['percent_trimmed'] = {
            'title': 'Trimmed',
            'description': '% Total Base Pairs trimmed',
            'max': 30,
            'min': 0,
            'scale': 'RdYlBu-rev',
            'format': '{:.1f}%'
        }
        self.general_stats_addcols(self.cutadapt_data, headers)
    

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
