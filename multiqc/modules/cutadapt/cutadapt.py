#!/usr/bin/env python

""" MultiQC module to parse output from Cutadapt """

from __future__ import print_function
from collections import OrderedDict
from distutils.version import StrictVersion
import logging
import re

from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Cutadapt module class, parses stderr logs.
    Also understands logs saved by Trim Galore!
    (which contain cutadapt logs)
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Cutadapt', anchor='cutadapt',
        href='https://cutadapt.readthedocs.io/',
        info="is a tool to find and remove adapter sequences, primers, poly-A"\
         "tails and other types of unwanted sequence from your high-throughput"\
         " sequencing reads.")

        # Find and load any Cutadapt reports
        self.cutadapt_data = dict()
        self.cutadapt_length_counts = dict()
        self.cutadapt_length_exp = dict()
        self.cutadapt_length_obsexp = dict()

        for f in self.find_log_files('cutadapt', filehandles=True):
            self.parse_cutadapt_logs(f)

        # Filter to strip out ignored sample names
        self.cutadapt_data = self.ignore_samples(self.cutadapt_data)

        if len(self.cutadapt_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.cutadapt_data)))

        # Write parsed report data to a file
        self.write_data_file(self.cutadapt_data, 'multiqc_cutadapt')

        # Basic Stats Table
        self.cutadapt_general_stats_table()

        # Bar plot with number of reads trimmed
        self.cutadapt_filtered_barplot()

        # Trimming Length Profiles
        self.cutadapt_length_trimmed_plot()


    def parse_cutadapt_logs(self, f):
        """ Go through log file looking for cutadapt output """
        fh = f['f']
        regexes = {
            '1.7': {
                'bp_processed': "Total basepairs processed:\s*([\d,]+) bp",
                'bp_written': "Total written \(filtered\):\s*([\d,]+) bp",
                'quality_trimmed': "Quality-trimmed:\s*([\d,]+) bp",
                'r_processed': "Total reads processed:\s*([\d,]+)",
                'pairs_processed': "Total read pairs processed:\s*([\d,]+)",
                'r_with_adapters': "Reads with adapters:\s*([\d,]+)",
                'r1_with_adapters': "Read 1 with adapter:\s*([\d,]+)",
                'r2_with_adapters': "Read 2 with adapter:\s*([\d,]+)",
                'r_too_short': "Reads that were too short:\s*([\d,]+)",
                'pairs_too_short': "Pairs that were too short:\s*([\d,]+)",
                'r_too_long': "Reads that were too long:\s*([\d,]+)",
                'pairs_too_long': "Pairs that were too long:\s*([\d,]+)",
                'r_too_many_N': "Reads with too many N:\s*([\d,]+)",
                'pairs_too_many_N': "Pairs with too many N:\s*([\d,]+)",
                'r_written': "Reads written \(passing filters\):\s*([\d,]+)",
                'pairs_written': "Pairs written \(passing filters\):\s*([\d,]+)",
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
        log_section = None
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
                # Manage case where sample name is '-' (reading from stdin)
                if s_name == '-':
                    s_name = f['s_name']
                else:
                    s_name = self.clean_s_name(s_name, f['root'])
                if s_name in self.cutadapt_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.cutadapt_data[s_name] = dict()

            if s_name is not None:
                self.add_data_source(f, s_name)

                # Search regexes for overview stats
                for k, r in regexes[cutadapt_version].items():
                    match = re.search(r, l)
                    if match:
                        self.cutadapt_data[s_name][k] = int(match.group(1).replace(',', ''))

                # Starting a new section
                if '===' in l:
                    log_section = l.strip().strip('=').strip()

                # Histogram showing lengths trimmed
                if 'length' in l and 'count' in l and 'expect' in l:
                    plot_sname = s_name
                    if log_section is not None:
                        plot_sname = '{} - {}'.format(s_name, log_section)
                    self.cutadapt_length_counts[plot_sname] = dict()
                    self.cutadapt_length_exp[plot_sname] = dict()
                    self.cutadapt_length_obsexp[plot_sname] = dict()
                    # Nested loop to read this section while the regex matches
                    for l in fh:
                        r_seqs = re.search("^(\d+)\s+(\d+)\s+([\d\.]+)", l)
                        if r_seqs:
                            a_len = int(r_seqs.group(1))
                            self.cutadapt_length_counts[plot_sname][a_len] = int(r_seqs.group(2))
                            self.cutadapt_length_exp[plot_sname][a_len] = float(r_seqs.group(3))
                            if float(r_seqs.group(3)) > 0:
                                self.cutadapt_length_obsexp[plot_sname][a_len] = float(r_seqs.group(2)) / float(r_seqs.group(3))
                            else:
                                # Cheating, I know. Infinity is difficult to plot.
                                self.cutadapt_length_obsexp[plot_sname][a_len] = float(r_seqs.group(2))
                        else:
                            break

        # Calculate a few extra numbers of our own
        for s_name, d in self.cutadapt_data.items():
            # Percent trimmed
            if 'bp_processed' in d and 'bp_written' in d:
                self.cutadapt_data[s_name]['percent_trimmed'] = (float(d['bp_processed'] - d['bp_written']) / d['bp_processed']) * 100
            elif 'bp_processed' in d and 'bp_trimmed' in d:
                self.cutadapt_data[s_name]['percent_trimmed'] = ((float(d.get('bp_trimmed', 0)) + float(d.get('quality_trimmed', 0))) / d['bp_processed']) * 100

    def cutadapt_general_stats_table(self):
        """ Take the parsed stats from the Cutadapt report and add it to the
        basic stats table at the top of the report """

        headers = {}
        headers['percent_trimmed'] = {
            'title': '% BP Trimmed',
            'description': '% Total Base Pairs trimmed',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlBu-rev'
        }
        self.general_stats_addcols(self.cutadapt_data, headers)

    def cutadapt_filtered_barplot(self):
        """ Bar plot showing proportion of reads trimmed """

        pconfig = {
            'id': 'cutadapt_filtered_reads_plot',
            'title': 'Cutadapt: Filtered Reads',
            'ylab': 'Counts'
        }

        # We just use all categories. If a report is generated with a mixture
        # of SE and PE data then this means quite a lot of categories.
        # Usually, only a single data type is used though - in that case
        # any categories with 0 across all samples will be ignored.
        cats = OrderedDict()
        cats['pairs_written'] = { 'name': 'Pairs passing filters' }
        cats['r_written'] = { 'name': 'Reads passing filters' }
        cats['pairs_too_short'] = { 'name': 'Pairs that were too short' }
        cats['r_too_short'] = { 'name': 'Reads that were too short' }
        cats['pairs_too_long'] = { 'name': 'Pairs that were too long' }
        cats['r_too_long'] = { 'name': 'Reads that were too long' }
        cats['pairs_too_many_N'] = { 'name': 'Pairs with too many N' }
        cats['r_too_many_N'] = { 'name': 'Reads with too many N' }

        self.add_section(
            name = 'Filtered Reads',
            anchor = 'cutadapt_filtered_reads',
            description = 'This plot shows the number of reads (SE) / pairs (PE) removed by Cutadapt.',
            plot = bargraph.plot(self.cutadapt_data, cats, pconfig)
        )

    def cutadapt_length_trimmed_plot (self):
        """ Generate the trimming length plot """

        pconfig = {
            'id': 'cutadapt_trimmed_sequences_plot',
            'title': 'Cutadapt: Lengths of Trimmed Sequences',
            'ylab': 'Counts',
            'xlab': 'Length Trimmed (bp)',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [
                {'name': 'Counts', 'ylab': 'Count'},
                {'name': 'Obs/Exp', 'ylab': 'Observed / Expected'}
            ]
        }

        self.add_section(
            name = 'Trimmed Sequence Lengths',
            anchor = 'cutadapt_trimmed_sequences',
            description = 'This plot shows the number of reads with certain lengths of adapter trimmed.',
            helptext = '''
            Obs/Exp shows the raw counts divided by the number expected due to sequencing errors.
            A defined peak may be related to adapter length.

            See the [cutadapt documentation](http://cutadapt.readthedocs.org/en/latest/guide.html#how-to-read-the-report)
            for more information on how these numbers are generated.
            ''',
            plot = linegraph.plot([self.cutadapt_length_counts, self.cutadapt_length_obsexp], pconfig)
        )
