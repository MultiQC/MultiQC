#!/usr/bin/env python

""" MultiQC module to parse output from Cutadapt """

from __future__ import print_function
import logging
import re

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

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
        
        for f in self.find_log_files(config.sp['cutadapt'], filehandles=True):
            self.parse_cutadapt_logs(f)        

        if len(self.cutadapt_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.cutadapt_data)))

        # Write parsed report data to a file
        self.write_data_file(self.cutadapt_data, 'multiqc_cutadapt')

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
            '1.7': {
                'bp_processed': "Total basepairs processed:\s*([\d,]+) bp",
                'bp_written': "Total written \(filtered\):\s*([\d,]+) bp",
                'quality_trimmed': "Quality-trimmed:\s*([\d,]+) bp",
                'r_processed': "Total reads processed:\s*([\d,]+)",
                'r_with_adapters': "Reads with adapters:\s*([\d,]+)"
            },
            '1.6': {
                'r_processed': "Processed reads:\s*([\d,]+)",
                'bp_processed': "Processed bases:\s*([\d,]+) bp",
                'r_trimmed': "Trimmed reads:\s*([\d,]+)",
                'bp_trimmed': "Trimmed bases:\s*([\d,]+) bp",
                'too_short': "Too short reads:\s*([\d,]+)",
                'too_long': "Too long reads:\s*([\d,]+)",
            }
        }
        s_name = None
        cutadapt_version = '1.7'
        for l in fh:
            # New log starting
            if l.startswith('This is cutadapt'):
                s_name = None
                c_version = re.match(r'This is cutadapt ([\d\.]+)', l)
                if c_version:
                    try:
                        assert(float(c_version.group(1)) <= 1.6)
                        cutadapt_version = '1.6'
                    except:
                        cutadapt_version = '1.7'
            
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
                self.add_data_source(f, s_name)
                
                # Search regexes for overview stats
                for k, r in regexes[cutadapt_version].items():
                    match = re.search(r, l)
                    if match:
                        self.cutadapt_data[s_name][k] = int(match.group(1).replace(',', ''))
                
                # Histogram showing lengths trimmed
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
        for s_name, d in self.cutadapt_data.items():
            if 'bp_processed' in d and 'bp_written' in d:
                self.cutadapt_data[s_name]['percent_trimmed'] = (float(d['bp_processed'] - d['bp_written']) / d['bp_processed']) * 100
            elif 'bp_processed' in d and 'bp_trimmed' in d:
                self.cutadapt_data[s_name]['percent_trimmed'] = (float(d['bp_trimmed']) / d['bp_processed']) * 100



    def cutadapt_general_stats_table(self):
        """ Take the parsed stats from the Cutadapt report and add it to the
        basic stats table at the top of the report """

        headers = {}
        headers['percent_trimmed'] = {
            'title': '% Trimmed',
            'description': '% Total Base Pairs trimmed',
            'max': 100,
            'min': 0,
            'suffix': '%',
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
            'title': 'Lengths of Trimmed Sequences',
            'ylab': 'Counts',
            'xlab': 'Length Trimmed (bp)',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [{'name': 'Counts', 'ylab': 'Count'},
                            {'name': 'Obs/Exp', 'ylab': 'Observed / Expected'}]
        }
        
        html += plots.linegraph.plot([self.cutadapt_length_counts, self.cutadapt_length_obsexp], pconfig)
        
        return html
