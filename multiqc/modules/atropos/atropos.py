#!/usr/bin/env python

""" MultiQC module to parse output from Atropos """

from __future__ import print_function
import logging
import re
from distutils.version import StrictVersion

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    atropos module class, parses stderr logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='atropos', anchor='atropos',
        href='https://atropos.readthedocs.io/',
        info="is a tool to find and remove adapter sequences, primers, poly-A"\
         "tails and other types of unwanted sequence from your high-throughput"\
         " sequencing reads.")

        # Find and load any atropos reports
        self.atropos_data = dict()
        self.atropos_length_counts = dict()
        self.atropos_length_exp = dict()
        self.atropos_length_obsexp = dict()

        for f in self.find_log_files('atropos', filehandles=True):
            self.parse_atropos_logs(f)

        # Filter to strip out ignored sample names
        self.atropos_data = self.ignore_samples(self.atropos_data)

        if len(self.atropos_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.atropos_data)))

        # Write parsed report data to a file
        self.write_data_file(self.atropos_data, 'multiqc_atropos')

        # Basic Stats Table
        self.atropos_general_stats_table()

        # Trimming Length Profiles
        self.atropos_length_trimmed_plot()


    def parse_atropos_logs(self, f):
        """ Go through log file looking for atropos output """
        fh = f['f']
        regexes = {
                'bp_processed': "Total bp processed:\s*([\d,]+) bp",
                'bp_written': "Total bp written:\s*([\d,]+) bp",
                'quality_trimmed': "Quality-trimmed:\s*([\d,]+) bp",
                'r_processed': "Total reads processed:\s*([\d,]+)",
                'r_written': "Reads written \(passing filters\):\s*([\d]+)",
                'r_with_adapters': "Reads with adapter:\s*([\d,]+)"
        }
        s_name = None
        for l in fh:
            # New log starting
            if 'Atropos' in l:
                s_name = None
                atropos_version = re.match(r'This is Atropos ([\d\.]+)', l)

            # Get sample name from end of command line params
            if l.startswith('Command line parameters'):
                s_name = l.split()[-1]
                # Manage case where sample name is '-' (reading from stdin)
                if s_name == '-':
                    s_name = f['s_name']
                else:
                    s_name = self.clean_s_name(s_name, f['root'])
                if s_name in self.atropos_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.atropos_data[s_name] = dict()
                self.atropos_length_counts[s_name] = dict()
                self.atropos_length_exp[s_name] = dict()
                self.atropos_length_obsexp[s_name] = dict()

            if s_name is not None:
                self.add_data_source(f, s_name)

                # Search regexes for overview stats
                for k, r in regexes.items():
                    match = re.search(r, l)
                    if match:
                        self.atropos_data[s_name][k] = int(match.group(1).replace(',', ''))

                # Histogram showing lengths trimmed
                if 'length' in l and 'count' in l and 'expect' in l:
                    # Nested loop to read this section while the regex matches
                    for l in fh:
                        r_seqs = re.search("\s+(\d+)\s+(\d+,\d+)\s+([\d\.]+)", l)
                        if r_seqs:
                            a_len = int(r_seqs.group(1))
                            self.atropos_length_counts[s_name][a_len] = int(r_seqs.group(2).replace(",", ""))
                            self.atropos_length_exp[s_name][a_len] = float(r_seqs.group(3))
                            if float(r_seqs.group(3)) > 0:
                                self.atropos_length_obsexp[s_name][a_len] = float(r_seqs.group(2).replace(",", "")) / float(r_seqs.group(3))
                            else:
                                # Cheating, I know. Infinity is difficult to plot.
                                self.atropos_length_obsexp[s_name][a_len] = float(r_seqs.group(2).replace(",", ""))

        # Calculate a few extra numbers of our own
        for s_name, d in self.atropos_data.items():
            if 'r_processed' in d and 'r_written' in d:
                self.atropos_data[s_name]['percent_trimmed'] = (float(d['r_processed'] - d['r_written']) / d['r_processed']) * 100



    def atropos_general_stats_table(self):
        """ Take the parsed stats from the atropos report and add it to the
        basic stats table at the top of the report """

        headers = {}
        headers['percent_trimmed'] = {
            'title': '% Trimmed',
            'description': '% Total Reads trimmed',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlBu-rev'
        }
        self.general_stats_addcols(self.atropos_data, headers)


    def atropos_length_trimmed_plot (self):
        """ Generate the trimming length plot """

        description = 'This plot shows the number of reads with certain lengths of adapter trimmed. \n\
        Obs/Exp shows the raw counts divided by the number expected due to sequencing errors. A defined peak \n\
        may be related to adapter length. See the \n\
        <a href="http://atropos.readthedocs.org/en/latest/guide.html#how-to-read-the-report" target="_blank">atropos documentation</a> \n\
        for more information on how these numbers are generated.'

        pconfig = {
            'id': 'atropos_plot',
            'title': 'Lengths of Trimmed Sequences',
            'ylab': 'Counts',
            'xlab': 'Length Trimmed (bp)',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [{'name': 'Counts', 'ylab': 'Count'},
                            {'name': 'Obs/Exp', 'ylab': 'Observed / Expected'}]
        }

        self.add_section(
            description = description,
            plot = linegraph.plot([self.atropos_length_counts, self.atropos_length_obsexp], pconfig)
        )
