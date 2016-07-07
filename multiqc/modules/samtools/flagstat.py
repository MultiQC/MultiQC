# coding: utf-8
#!/usr/bin/env python

""" MultiQC submodule to parse output from Samtools flagstat """

import logging
from collections import OrderedDict
from multiqc import config, plots 

# Initialise the logger
log = logging.getLogger(__name__)

# flagstat has one thing per line, documented here (search for flagstat):
# http://www.htslib.org/doc/samtools.html 
# The numbers are line numbers for the file
LABELS= {1: 'total', 
        2: 'secondary',
        3: 'supplementary', 
        4: 'duplicates', 
        5: 'mapped', 
        6: 'paired in sequencing', 
        7: 'read1', 
        8: 'read2', 
        9: 'properly paired', 
        10: 'with itself and mate mapped', 
        11: 'singletons',
        12: 'with mate mapped to a different chr', 
        13: 'with make mapped to a different chr (mapQ >= 5)', 
        }

""" Take a filename, parse the data assuming it's a flagstat file
    Returns a dictionary {'lineName' : {'pass':value, 'fail':value}"""
def parse_single_report(file_thing):
    parsed_data = {}
    lines = file_thing.splitlines()
    for i, line in enumerate(lines, 1):
        d = {} # contains {'pass' : first number, 'fail' : second number}
        words = line.split(' ')
        d['pass'] = int(words[0].strip())
        d['fail'] = int(words[2].strip())
        parsed_data[LABELS[i]] = d
    return parsed_data

def parse_reports(self):
    """ Find Samtools flagstat logs and parse their data """

    self.samtools_flagstat = dict()
    for f in self.find_log_files(config.sp['samtools']['flagstat']):
        parsed_data = parse_single_report(f['f'])
        if len(parsed_data) > 0:
            if f['s_name'] in self.samtools_flagstat:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f)
            self.samtools_flagstat[f['s_name']] = parsed_data

    if len(self.samtools_flagstat) > 0:

        # Write parsed report data to a file
        print "lalala"
        self.write_data_file(self.samtools_flagstat, 'multiqc_samtools_flagstat')

        # General Stats Table
        self.general_flagstat_headers['error_rate'] = {
            'title': 'Error rate',
            'description': 'Error rate using CIGAR',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd',
            'format': '{:.2f}%',
            'modify': lambda x: x * 100.0
        }
        self.general_flagstat_headers['non-primary_alignments'] = {
            'title': 'M Non-Primary',
            'description': 'Non-primary alignments (millions)',
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        self.general_flagstat_headers['reads_mapped'] = {
            'title': 'M Reads Mapped',
            'description': 'Reads Mapped in the bam file',
            'min': 0,
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        self.general_flagstat_headers['reads_mapped_percent'] = {
            'title': '% Mapped',
            'description': '% Mapped Reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn',
            'format': '{:.1f}%'
        }
        self.general_flagstat_headers['raw_total_sequences'] = {
            'title': 'M Total seqs',
            'description': 'Total sequences in the bam file',
            'min': 0,
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        for s_name in self.samtools_flagstat:
            if s_name not in self.general_flagstat_data:
                self.general_flagstat_data[s_name] = dict()
            self.general_flagstat_data[s_name].update( self.samtools_flagstat[s_name] )
        
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
        keys['raw_total_sequences'] = dict(reads, **{'title': 'Total sequences' })
        keys['reads_mapped'] = dict(reads, **{'title': 'Mapped reads' })
        keys['reads_mapped_and_paired'] = dict(reads, **{'title': 'Mapped &amp; paired', 'description': 'Paired-end technology bit set + both mates mapped' })
        keys['reads_properly_paired'] = dict(reads, **{'title': 'Properly paired', 'description': 'Proper-pair bit set' })
        keys['reads_duplicated'] = dict(reads, **{'title': 'Duplicated', 'description': 'PCR or optical duplicate bit set' })
        keys['reads_unmapped'] = dict(reads, **{'title': 'Unmapped reads' })
        keys['reads_QC_failed'] = dict(reads, **{'title': 'QC Failed'})
        keys['reads_MQ0'] = dict(reads, **{'title': 'Reads MQ0', 'description': 'Reads mapped and MQ=0' })
        keys['bases_mapped_(cigar)'] = dict(bases, **{'title': 'Mapped bases (cigar)', 'description': 'Mapped bases (cigar)' })
        keys['bases_trimmed'] = dict(bases, **{'title': 'Bases Trimmed' })
        keys['bases_duplicated'] = dict(bases, **{'title': 'Duplicated bases' })
        keys['pairs_on_different_chromosomes'] = dict(reads, **{'title': 'Diff chromosomes', 'description': 'Pairs on different chromosomes' })
        keys['pairs_with_other_orientation'] = dict(reads, **{'title': 'Other orientation', 'description': 'Pairs with other orientation' })
        keys['inward_oriented_pairs'] = dict(reads, **{'title': 'Inward pairs', 'description': 'Inward oriented pairs' })
        keys['outward_oriented_pairs'] = dict(reads, **{'title': 'Outward pairs', 'description': 'Outward oriented pairs' })
        
        self.sections.append({
            'name': 'Samtools flagstat output',
            'anchor': 'samtools-flagstat',
            'content': '<p>This module parses the output from <code>samtools flagstat</code>. All numbers in millions.</p>' +
                        plots.beeswarm.plot(self.samtools_flagstat, keys, {'id': 'samtools-flagstat-dp'})
        })
    
    # Return the number of logs that were found
    return len(self.samtools_flagstat)


