#!/usr/bin/env python

""" MultiQC submodule to parse output from Rockhopper summary files
https://cs.wellesley.edu/~btjaden/Rockhopper/ """

from __future__ import print_function
import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialize the parent object
        super(MultiqcModule, self).__init__(
            name = 'Rockhopper',
            anchor = 'rockhopper',
            href = 'https://cs.wellesley.edu/~btjaden/Rockhopper/',
            info = """
            is a comprehensive and user-friendly system
            for computational analysis of bacterial RNA-seq data.
            It can align reads to genomes, assemble transcripts,
            identify transcript boundaries, and discover novel
            transcripts such as small RNAs.
            """
        )

        # Set up vars
        self.rh_data = dict()

        # Parse summary file
        for f in self.find_log_files('rockhopper'):
            self.parse_rockhopper_summary(f)

        # Filter to strip out ignored sample names
        self.rh_data = self.ignore_samples(self.rh_data)

        if len(self.rh_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.rh_data)))

        # Write to file
        self.write_data_file(self.rh_data, 'multiqc_rockhopper')

        # Basic stats Table
        self.rockhopper_general_stats_table()

        # Alignment bar plot
        self.rockhopper_count_bar_plot()


    def parse_rockhopper_summary(self,f):

        s_name = None

        # Initialize stats fields
        stats_index = [
            'mRNA-sense',
            'mRNA-antisense',
            'rRNA-sense',
            'rRNA-antisense',
            'tRNA-sense',
            'tRNA-antisense',
            'ncRNA-sense',
            'ncRNA-antisense',
            'unannotated'
        ]

        results = { name:0 for name in stats_index }

        # Parse rockhopper output line-by-line since there may be many genomes
        lines = f['f'].split('\n')
        for i, line in enumerate(lines):

            # Get sample name
            if line.startswith('Aligning sequencing reads from file:'):
                s_name = line.split(':', 1)[1].strip()
                s_name = self.clean_s_name(s_name,f['root'])

            elif line.startswith('Aligning sequencing reads from files:'):
                s_name = lines[i+1].strip()
                s_name = self.clean_s_name(s_name,f['root'])

            # Get total number of reads read by rockhopper
            if line.startswith('Total reads:'):
                results['total-reads'] = int(re.search('Total reads:\s*(\d*)',line).group(1))
                i += 1

            # Get number of reads aligned to each genome
            elif line.startswith('Successfully aligned reads'):

                # Get Genome ID
                genome = re.search('\(>(.*?) .*\)',line).group(1)

                # Get number of aligned reads
                genome_reads = int(re.search('Successfully aligned reads:\s*(\d*)',line).group(1))

                # Get percent of reads in each category
                stats = [int(re.search('(\d+)\%',subline).group(1)) for subline in lines[i+1:i+10]]
                for name,val in zip(stats_index,stats):
                    # Convert percentages to true number of reads in each category
                    results[name] += int(round(val*genome_reads/100))

                # Skip 10 lines
                i += 10

            else:
                i += 1

        if len(results) > 0 and s_name is not None:

            # Calculate unaligned read count
            total_mapped_reads  = sum([v for k,v in results.items() if k != 'total-reads'])
            results['unaligned'] = results['total-reads'] - total_mapped_reads

            if s_name in self.rh_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.rh_data[s_name] = results

    def rockhopper_general_stats_table(self):
        """ Take the parsed stats from the Rockhopper summary and add it to the
        basic stats table at the top of the report """

        headers = {}
        headers['mRNA-sense'] = {
            'title': 'CDS Reads ({})'.format(config.read_count_prefix),
            'description':'Reads aligned to coding regions ({})'.format(config.read_count_desc),
            'min':0,
            'scale':'Blues',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        headers['mRNA-antisense'] = {
            'title': 'CDS Reads (a/s, {})'.format(config.read_count_prefix),
            'description':'Antisense reads aligned to coding regions ({})'.format(config.read_count_desc),
            'min':0,
            'scale':'Blues',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count',
            'hidden':True
        }
        headers['rRNA-sense'] = {
            'title': 'rRNA Reads ({})'.format(config.read_count_prefix),
            'description':'Reads aligned to rRNA ({})'.format(config.read_count_desc),
            'min':0,
            'scale':'Blues',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        headers['rRNA-antisense'] = {
            'title': 'rRNA Reads (a/s, {})'.format(config.read_count_prefix),
            'description':'Antisense reads aligned to rRNA ({})'.format(config.read_count_desc),
            'min':0,
            'scale':'Blues',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count',
            'hidden':True
        }
        self.general_stats_addcols(self.rh_data,headers)

    def rockhopper_count_bar_plot(self):
        """ Stacked bar plot showing counts of reads """

        pconfig = {
            'id':'rockhopper_reads_counts_plot',
            'title': 'Rockhopper: Alignment types',
            'ylab': 'Number of reads',
            'tt_percentage': False
        }

        # Plot bar graph of groups
        keys = OrderedDict()
        keys['mRNA-sense']      = {'name': "mRNA (Sense)"}
        keys['mRNA-antisense']  = {'name': "mRNA (Antisense)"}
        keys['rRNA-sense']      = {'name': "rRNA (Sense)"}
        keys['rRNA-antisense']  = {'name': "rRNA (Antisense)"}
        keys['tRNA-sense']      = {'name': "tRNA (Sense)"}
        keys['tRNA-antisense']  = {'name': "tRNA (Antisense)"}
        keys['ncRNA-sense']     = {'name': "ncRNA (Sense)"}
        keys['ncRNA-antisense'] = {'name': "ncRNA (Antisense)"}
        keys['unannotated']     = {'name': "Unannotated"}
        keys['unaligned']       = {'name':'Unaligned'}

        self.add_section (
            name = 'Rockhopper',
            anchor = 'rockhopper',
            description = """
            This plot shows the number of reads mapped to
            different regions of the genome, accounting for
            sense and antisense alignment, if relevant.
            """,
            plot = bargraph.plot(self.rh_data, keys, pconfig)
        )
