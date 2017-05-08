#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from Samtools stats """

import logging
from collections import OrderedDict
from multiqc import config
from multiqc.plots import beeswarm, bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class StatsReportMixin():
    """ Mixin class, loaded by main samtools MuliqcModule class. """

    def parse_samtools_stats(self):
        """ Find Samtools stats logs and parse their data """

        self.samtools_stats = dict()
        for f in self.find_log_files('samtools/stats'):
            parsed_data = dict()
            for line in f['f'].splitlines():
                if not line.startswith("SN"):
                    continue
                sections = line.split("\t")
                field = sections[1].strip()[:-1]
                field = field.replace(' ', '_')
                value = float(sections[2].strip())
                parsed_data[field] = value

            if len(parsed_data) > 0:
                # Work out some percentages
                if 'raw_total_sequences' in parsed_data:
                    for k in list(parsed_data.keys()):
                        if k.startswith('reads_') and k != 'raw_total_sequences' and parsed_data['raw_total_sequences'] > 0:
                            parsed_data['{}_percent'.format(k)] = (parsed_data[k] / parsed_data['raw_total_sequences']) * 100

                if f['s_name'] in self.samtools_stats:
                    log.debug("Duplicate sample name found! Overwriting: {}"
                              .format(f['s_name']))
                self.add_data_source(f, section='stats')
                self.samtools_stats[f['s_name']] = parsed_data

        # Filter to strip out ignored sample names
        self.samtools_stats = self.ignore_samples(self.samtools_stats)

        if len(self.samtools_stats) > 0:

            # Write parsed report data to a file
            self.write_data_file(self.samtools_stats, 'multiqc_samtools_stats')

            # General Stats Table
            stats_headers = OrderedDict()
            stats_headers['error_rate'] = {
                'title': 'Error rate',
                'description': 'Error rate using CIGAR',
                'min': 0,
                'max': 100,
                'suffix': '%',
                'scale': 'OrRd',
                'format': '{:,.2f}',
                'modify': lambda x: x * 100.0
            }
            stats_headers['non-primary_alignments'] = {
                'title': '{} Non-Primary'.format(config.read_count_prefix),
                'description': 'Non-primary alignments ({})'.format(config.read_count_desc),
                'min': 0,
                'scale': 'PuBu',
                'modify': lambda x: x * config.read_count_multiplier,
                'shared_key': 'read_count'
            }
            stats_headers['reads_mapped'] = {
                'title': '{} Reads Mapped'.format(config.read_count_prefix),
                'description': 'Reads Mapped in the bam file ({})'.format(config.read_count_desc),
                'min': 0,
                'modify': lambda x: x * config.read_count_multiplier,
                'shared_key': 'read_count'
            }
            stats_headers['reads_mapped_percent'] = {
                'title': '% Mapped',
                'description': '% Mapped Reads',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'RdYlGn'
            }
            stats_headers['raw_total_sequences'] = {
                'title': '{} Total seqs'.format(config.read_count_prefix),
                'description': 'Total sequences in the bam file ({})'.format(config.read_count_desc),
                'min': 0,
                'modify': lambda x: x * config.read_count_multiplier,
                'shared_key': 'read_count'
            }
            self.general_stats_addcols(self.samtools_stats, stats_headers, 'Samtools Stats')

            # Make bargraph plot of mapped/unmapped reads
            self.alignment_section(self.samtools_stats)

            # Make dot plot of counts
            keys = OrderedDict()
            reads = {
                'min': 0,
                'modify': lambda x: float(x) / 1000000.0,
                'suffix': 'M reads',
                'decimalPlaces': 2,
                'shared_key': 'read_count'
            }
            bases = {
                'min': 0,
                'modify': lambda x: float(x) / 1000000.0,
                'suffix': 'M bases',
                'decimalPlaces': 2,
                'shared_key': 'base_count'
            }
            keys['raw_total_sequences'] = dict(reads, **{'title': 'Total sequences'})
            keys['reads_mapped_and_paired'] = dict(reads, **{'title': 'Mapped &amp; paired', 'description': 'Paired-end technology bit set + both mates mapped' })
            keys['reads_properly_paired'] = dict(reads, **{'title': 'Properly paired', 'description': 'Proper-pair bit set'})
            keys['reads_duplicated'] = dict(reads, **{'title': 'Duplicated', 'description': 'PCR or optical duplicate bit set'})
            keys['reads_QC_failed'] = dict(reads, **{'title': 'QC Failed'})
            keys['reads_MQ0'] = dict(reads, **{'title': 'Reads MQ0', 'description': 'Reads mapped and MQ=0'})
            keys['bases_mapped_(cigar)'] = dict(bases, **{'title': 'Mapped bases (cigar)', 'description': 'Mapped bases (cigar)'})
            keys['bases_trimmed'] = dict(bases, **{'title': 'Bases Trimmed'})
            keys['bases_duplicated'] = dict(bases, **{'title': 'Duplicated bases'})
            keys['pairs_on_different_chromosomes'] = dict(reads, **{'title': 'Diff chromosomes', 'description': 'Pairs on different chromosomes'})
            keys['pairs_with_other_orientation'] = dict(reads, **{'title': 'Other orientation', 'description': 'Pairs with other orientation'})
            keys['inward_oriented_pairs'] = dict(reads, **{'title': 'Inward pairs', 'description': 'Inward oriented pairs'})
            keys['outward_oriented_pairs'] = dict(reads, **{'title': 'Outward pairs', 'description': 'Outward oriented pairs'})

            self.add_section (
                name = 'Alignment metrics',
                anchor = 'samtools-stats',
                description = "This module parses the output from <code>samtools stats</code>. All numbers in millions.",
                plot = beeswarm.plot(self.samtools_stats, keys, {'id': 'samtools-stats-dp'})
            )

        # Return the number of logs that were found
        return len(self.samtools_stats)


    def alignment_section(self, samples_data):
        bedgraph_data = {}
        for sample_id, data in samples_data.items():
            expected_total = data['raw_total_sequences']
            read_sum = (data['reads_mapped'] + data['reads_unmapped'])
            if read_sum == expected_total:
                bedgraph_data[sample_id] = data
            else:
                log.warn("sum of mapped/unmapped reads not matching total, "
                         "skipping samtools plot for: {}".format(sample_id))
        self.add_section (
            name = 'Percent Mapped',
            anchor = 'samtools-stats-alignment',
            description = "Alignment metrics from <code>samtools stats</code>; mapped vs. unmapped reads.",
            plot = alignment_chart(bedgraph_data)
        )


def alignment_chart(data):
    """Make the HighCharts HTML to plot the alignment rates """
    keys = OrderedDict()
    keys['reads_mapped'] = {'color': '#437bb1', 'name': 'Mapped'}
    keys['reads_unmapped'] = {'color': '#b1084c', 'name': 'Unmapped'}

    # Config for the plot
    plot_conf = {
        'id': 'samtools_alignment_plot',
        'title': 'Samtools Alignment Scores',
        'ylab': '# Reads',
        'cpswitch_counts_label': 'Number of Reads'
    }
    return bargraph.plot(data, keys, plot_conf)
