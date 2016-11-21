#!/usr/bin/env python

""" MultiQC module to parse output from Slamdunk """

from __future__ import print_function
import logging
import re
from distutils.version import StrictVersion
from collections import OrderedDict

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Slamdunk module class, parses slamdunk logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Slamdunk', anchor='slamdunk',
        href='https://github.com/t-neumann/slamdunk',
        info="is a tool to analyze SLAMSeq data."\
         "blabla"\
         " blabla.")

        self.slamdunk_data = dict()
        self.rates_data = dict()
#         # Find and load any Cutadapt reports
#         self.cutadapt_data = dict()
#         self.cutadapt_length_counts = dict()
#         self.cutadapt_length_exp = dict()
#         self.cutadapt_length_obsexp = dict()
        
        
        for f in self.find_log_files(config.sp['slamdunk']['utrrates']):
            #self.parse_cutadapt_logs(f)
            #log.info(f['f'])
            log.info(f['s_name'])
            log.info(f['root'])
            log.info(f['fn'])
            
        for f in self.find_log_files(config.sp['slamdunk']['summary'], filehandles = True):
            self.parseSummary(f)
            
        for f in self.find_log_files(config.sp['slamdunk']['rates'], filehandles = True):
            self.parseSlamdunkRates(f)

        if len(self.slamdunk_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.slamdunk_data)))

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.SlamdunkGeneralStatsTable()

        # Trimming Length Profiles
        # Only one section, so add to the intro
#         self.intro += self.cutadapt_length_trimmed_plot()

    def parseSlamdunkRates(self, f):
        for line in f['f']:
            log.info(line)
#                 if (not line.startswith('#') and not line.startswith('FileName')):
#                     fields = line.rstrip().split("\t")
#                     self.slamdunk_data[self.clean_s_name(fields[0],"")] = dict()
#                     self.slamdunk_data[self.clean_s_name(fields[0],"")]['sequenced'] = fields[4]
#                     self.slamdunk_data[self.clean_s_name(fields[0],"")]['mapped'] = fields[5]
#                     self.slamdunk_data[self.clean_s_name(fields[0],"")]['deduplicated'] = fields[6]
#                     self.slamdunk_data[self.clean_s_name(fields[0],"")]['filtered'] = fields[7]

    
    def parseSummary(self, f):
        for line in f['f']:
            if (not line.startswith('#') and not line.startswith('FileName')):
                fields = line.rstrip().split("\t")
                self.slamdunk_data[self.clean_s_name(fields[0],"")] = dict()
                self.slamdunk_data[self.clean_s_name(fields[0],"")]['sequenced'] = int(fields[4])
                self.slamdunk_data[self.clean_s_name(fields[0],"")]['mapped'] = int(fields[5])
                self.slamdunk_data[self.clean_s_name(fields[0],"")]['deduplicated'] = int(fields[6])
                self.slamdunk_data[self.clean_s_name(fields[0],"")]['filtered'] = int(fields[7])
        self.add_data_source(f)


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
                'quality_trimmed': "Quality-trimmed:\s*([\d,]+) bp",
                'bp_trimmed': "Trimmed bases:\s*([\d,]+) bp",
                'too_short': "Too short reads:\s*([\d,]+)",
                'too_long': "Too long reads:\s*([\d,]+)",
            }
        }
        s_name = None
        cutadapt_version = '1.7'
        for l in fh:
            # New log starting
            if 'cutadapt' in l:
                s_name = None
                c_version = re.match(r'This is cutadapt ([\d\.]+)', l)
                if c_version:
                    try:
                        assert(StrictVersion(c_version.group(1)) <= StrictVersion('1.6'))
                        cutadapt_version = '1.6'
                    except:
                        cutadapt_version = '1.7'
                c_version_old = re.match(r'cutadapt version ([\d\.]+)', l)
                if c_version_old:
                    try:
                        assert(StrictVersion(c_version.group(1)) <= StrictVersion('1.6'))
                        cutadapt_version = '1.6'
                    except:
                        # I think the pattern "cutadapt version XX" is only pre-1.6?
                        cutadapt_version = '1.6'
            
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
                self.cutadapt_data[s_name]['percent_trimmed'] = ((float(d.get('bp_trimmed', 0)) + float(d.get('quality_trimmed', 0))) / d['bp_processed']) * 100



    def SlamdunkGeneralStatsTable(self):
        """ Take the parsed summary stats from Slamdunk and add it to the
        basic stats table at the top of the report """


        headers = OrderedDict()
        headers['sequenced'] = {
            'title': 'Sequenced',
            'description': '# sequenced reads',
            'shared_key': 'reads',
            'min': 0,
            'format': '{:.f}'
        }
        headers['mapped'] = {
            'title': 'Mapped',
            'description': '# mapped reads',
            'shared_key': 'reads',
            'min': 0,
            'format': '{:.f}'
        }
        headers['deduplicated'] = {
            'title': 'Deduplicated',
            'description': '# deduplicated reads',
            'shared_key': 'reads',
            'min': 0,
            'format': '{:.f}'
        }
        headers['filtered'] = {
            'title': 'Filtered',
            'description': '# reads after filtering',
            'shared_key': 'reads',
            'min': 0,
            'format': '{:.f}'
        }
        self.general_stats_addcols(self.slamdunk_data, headers)
    

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
