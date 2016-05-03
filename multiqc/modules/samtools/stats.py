#!/usr/bin/env python

""" MultiQC submodule to parse output from Samtools stats """

import logging
from collections import OrderedDict
from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Samtools stats logs and parse their data """

    self.samtools_stats = dict()
    for f in self.find_log_files(config.sp['samtools']['stats']):
        parsed_data = dict()
        for line in f['f'].splitlines():
            if not line.startswith("SN"):
                continue
            sections = line.split("\t")
            field = sections[1].strip()[:-1]
            field = field.replace(' ','_')
            value = float(sections[2].strip())
            parsed_data[field] = value
        if len(parsed_data) > 0:
            if f['s_name'] in self.samtools_stats:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f)
            self.samtools_stats[f['s_name']] = parsed_data

    if len(self.samtools_stats) > 0:

        # Write parsed report data to a file
        self.write_data_file(self.samtools_stats, 'multiqc_samtools_stats')

        # General Stats Table
        self.general_stats_headers['raw_total_sequences'] = {
            'title': 'M Total seqs',
            'description': 'Total sequences in the bam file',
            'min': 0,
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_headers['reads_mapped'] = {
            'title': 'M Reads Mapped',
            'description': 'Reads Mapped in the bam file',
            'min': 0,
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_headers['non-primary_alignments'] = {
            'title': 'M Non-Primary Alignments',
            'description': 'Non primary alignment (millions)',
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_headers['error_rate'] = {
            'title': 'Error rate',
            'description': 'Error rate using CIGAR',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd',
            'format': '{:.2f}%',
            'modify': lambda x: x * 100.0
        }
        for s_name in self.samtools_stats:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.samtools_stats[s_name] )
        
        # Make dot plot of counts
        keys = OrderedDict()
        reads = {
            'min': 0,
            'modify': lambda x: float(x) / 1000000.0,
            'decimalPlaces': 2,
            'shared_key': 'read_count'
        }
        bases = {
            'min': 0,
            'modify': lambda x: float(x) / 1000000.0,
            'decimalPlaces': 2,
            'shared_key': 'base_count'
        }
        keys['raw_total_sequences'] = dict(reads, **{'title': 'Total sequences'})
        keys['reads_mapped'] = dict(reads, **{'title': 'Mapped reads'})
        keys['reads_mapped_and_paired'] = dict(reads, **{'title': 'Mapped &amp; paired'})
        keys['reads_properly_paired'] = dict(reads, **{'title': 'Properly paired'})
        keys['reads_duplicated'] = dict(reads, **{'title': 'Duplicated'})
        keys['reads_unmapped'] = dict(reads, **{'title': 'Unmapped reads'})
        keys['reads_QC_failed'] = dict(reads, **{'title': 'QC Failed'})
        keys['reads_MQ0'] = dict(reads, **{'title': 'Reads MQ0'})
        keys['bases_mapped'] = dict(bases, **{'title': 'Mapped bases'})
        keys['bases_mapped_(cigar)'] = dict(bases, **{'title': 'Mapped bases (cigar)'})
        keys['bases_trimmed'] = dict(bases, **{'title': 'Bases Trimmed'})
        keys['bases_duplicated'] = dict(bases, **{'title': 'Duplicated bases'})
        keys['pairs_on_different_chromosomes'] = dict(reads, **{'title': 'Diff chromosomes', 'description': 'Pairs on different chromosomes'})
        keys['pairs_with_other_orientation'] = dict(reads, **{'title': 'Other orientation', 'description': 'Pairs with other orientation'})
        keys['inward_oriented_pairs'] = dict(reads, **{'title': 'Inward pairs', 'description': 'Inward oriented pairs' })
        keys['outward_oriented_pairs'] = dict(reads, **{'title': 'Outward pairs', 'description': 'Outward oriented pairs' })
        
        self.sections.append({
            'name': 'samtools Stats',
            'anchor': 'samtools-stats',
            'content': plots.beeswarm.plot(self.samtools_stats, keys)
        })
    
    # Return the number of logs that were found
    return len(self.samtools_stats)


