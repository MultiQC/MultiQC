#!/usr/bin/env python

""" MultiQC module to parse output from odgi stats """

from __future__ import print_function
import logging
from collections import OrderedDict

from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Some class desc #TODO
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Odgi', anchor='odgi',
                                            href='Some ODGI URL',  # TODO
                                            info="Some info string!")  # TODO

        self.odgi_stats_map = dict()

        # Find and load any odgi stats file
        log_files = self.find_log_files('odgi/stats')
        if not log_files:
            raise UserWarning
        for f in log_files:
            self.parse_odgi_stats_report(f)
            self.add_data_source(f)

        log.info("Found {} reports".format(len(self.odgi_stats_map)))

        # Write parsed report data to a file
        self.write_data_file(self.odgi_stats_map, 'multiqc_odgi_stats')

        self.plot_odgi_stats()

        # Assignment bar plot
        # self.rsem_mapped_reads_plot()

        # Multimapping line plot
        # self.rsem_multimapping_plot()

    def parse_odgi_stats_report(self, f):
        """
        Bla
        """
        lines = f['f'].splitlines()
        # parse the generals stats
        stats = lines[1].split('\t')
        length, nodes, edges, paths = stats[0], stats[1], stats[2], stats[3]
        # parse the mean links length
        mean_links_length = lines[4].split('\t')
        path, in_node_space, in_nucleotide_space, num_links_considered = mean_links_length[0], mean_links_length[1], mean_links_length[2], mean_links_length[3]
        # parse sum of path nodes distances
        sum_of_path_nodes_distances = lines[7].split('\t')
        path, in_node_space, in_nucleotide_space, \
        nodes, nucleotides, num_penalties, \
        num_penalties_different_orientation = sum_of_path_nodes_distances[0], sum_of_path_nodes_distances[1], \
                                              sum_of_path_nodes_distances[2], sum_of_path_nodes_distances[3], \
                                              sum_of_path_nodes_distances[4], sum_of_path_nodes_distances[5], \
                                              sum_of_path_nodes_distances[6]
        self.odgi_stats_map.update(MultiqcModule.compress_stats_data(f['fn'], stats, mean_links_length, sum_of_path_nodes_distances))

    def plot_odgi_stats(self):
        """
        Plot the odgi stats
        """
        data = dict()

        headers = OrderedDict()
        headers['Nodes'] = {
            'title': 'Nodes',
            'description': 'My Second Column',
            'max': 100,
            'min': 0,
            'scale': 'Set1',
        }
        headers['Edges'] = {
            'title': 'Edges',
            'description': 'My Second Column',
            'max': 100,
            'min': 0,
            'scale': 'Set1'
        }
        headers['Length'] = {
            'title': 'Length',
            'description': 'My Second Column',
            'max': 100,
            'min': 0,
            'scale': 'Set1'
        }
        headers['Paths'] = {
            'title': 'Paths',
            'description': 'My Second Column',
            'max': 100,
            'min': 0,
            'scale': 'Set1'
        }
        for i, fn in enumerate(self.odgi_stats_map.keys()):
            file_stats = self.odgi_stats_map[fn]
            data.update({f'sample{i}': file_stats[f'General stats {fn}']})
        log.info(data)
        self.general_stats_addcols(data, headers)

    @staticmethod
    def compress_stats_data(fn, stats, mean_links_length, sum_of_path_nodes_distances) -> dict:
        """
        Bla
        """
        return {
            fn: {
                f'General stats {fn}': {
                    'Length': stats[0],
                    'Nodes': stats[1],
                    'Edges': stats[2],
                    'Paths': stats[3],
                },
                'Mean_links_length': {
                    'Path': mean_links_length[0],
                    'In_node_space': mean_links_length[1],
                    'In_nucleotide_space': mean_links_length[2],
                    'Num_links_considered': mean_links_length[3]
                },
                'Sum_of_path_nodes_distances': {
                    'Path': sum_of_path_nodes_distances[0],
                    'In_node_space': sum_of_path_nodes_distances[1],
                    'In_nucleotide_space': sum_of_path_nodes_distances[2],
                    'Nodes': sum_of_path_nodes_distances[3],
                    'Nucleotides': sum_of_path_nodes_distances[4],
                    'Num_penalties': sum_of_path_nodes_distances[5],
                    'Num_penalties_different_orientation': sum_of_path_nodes_distances[6]
                }
            }
        }

    def rsem_mapped_reads_plot(self):
        """ Make the rsem assignment rates plot """

        # Plot categories
        keys = OrderedDict()
        keys['Unique'] = {'color': '#437bb1', 'name': 'Aligned uniquely to a gene'}
        keys['Multi'] = {'color': '#e63491', 'name': 'Aligned to multiple genes'}
        keys['Filtered'] = {'color': '#b1084c', 'name': 'Filtered due to too many alignments'}
        keys['Unalignable'] = {'color': '#7f0000', 'name': 'Unalignable reads'}

        # Config for the plot
        config = {
            'id': 'rsem_assignment_plot',
            'title': 'RSEM: Mapped reads',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads',
            'hide_zero_cats': False
        }

        self.add_section(
            name='Mapped Reads',
            anchor='rsem_mapped_reads',
            description='A breakdown of how all reads were aligned for each sample.',
            plot=bargraph.plot(self.rsem_mapped_data, keys, config)
        )

    def rsem_multimapping_plot(self):
        """ Make a line plot showing the multimapping levels """

        pconfig = {
            'id': 'rsem_multimapping_rates',
            'title': 'RSEM: Multimapping Rates',
            'ylab': 'Counts',
            'xlab': 'Number of alignments',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>{point.x} alignments</b>: {point.y:.0f}',
        }

        self.add_section(
            name='Multimapping rates',
            anchor='rsem_multimapping',
            description='A frequency histogram showing how many reads were aligned to `n` reference regions.',
            helptext='''In an ideal world, every sequence reads would align uniquely to a single location in the
                reference. However, due to factors such as repeititve sequences, short reads and sequencing errors,
                reads can be align to the reference 0, 1 or more times. This plot shows the frequency of each factor
                of multimapping. Good samples should have the majority of reads aligning once.''',
            plot=linegraph.plot(self.rsem_multimapping_data, pconfig)
        )
