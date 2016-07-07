# coding: utf-8
#!/usr/bin/env python

""" MultiQC submodule to parse output from Samtools flagstat """

import logging
from collections import OrderedDict, defaultdict
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
        13: 'with mate mapped to a different chr (mapQ >= 5)', 
        }

""" Take a filename, parse the data assuming it's a flagstat file
    Returns a dictionary {'lineName_pass' : value, 'lineName_fail' : value}"""
def parse_single_report(file_thing):
    parsed_data = {}
    lines = file_thing.splitlines()
    for i, line in enumerate(lines, 1):
        prefix = LABELS[i] + '_'
        words = line.split(' ')
        parsed_data[prefix + 'pass'] = int(words[0].strip())
        parsed_data[prefix + 'fail'] = int(words[2].strip())
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

        # Write parsed report data to a file (restructure first)
        to_write = defaultdict(dict)
        self.write_data_file(self.samtools_flagstat, 'multiqc_samtools_flagstat')

        """
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
        """
        
        # Make dot plot of counts
        keys = OrderedDict()
        reads = {
            'min': 0,
            'modify': lambda x: float(x) / 1000000.0,
            'suffix': 'M reads',
            'decimalPlaces': 2,
            'shared_key': 'read_count'
        }
        keys['total_pass'] = dict(reads, **{'title': 'Total Passed' })
        keys['total_fail'] = dict(reads, **{'title': 'Total Failed' })
        keys['secondary_pass'] = dict(reads, **{'title': 'Secondary Passed' })
        keys['secondary_fail'] = dict(reads, **{'title': 'Secondary Failed' })
        keys['supplementary_pass'] = dict(reads, **{'title': 'Supplementary Passed' })
        keys['supplementary_fail'] = dict(reads, **{'title': 'Supplementary Failed' })
        keys['duplicates_pass'] = dict(reads, **{'title': 'Duplicate Passed' })
        keys['duplicates_fail'] = dict(reads, **{'title': 'Duplicate Failed' })
        keys['mapped_pass'] = dict(reads, **{'title': 'Mapped Passed' })
        keys['mapped_fail'] = dict(reads, **{'title': 'Mapped Failed' })
        keys['paired in sequencing_pass'] = dict(reads, **{'title': 'Paired in Sequencing Passed' })
        keys['paired in sequencing_fail'] = dict(reads, **{'title': 'Paired in Sequencing Failed' })
        keys['read1_pass'] = dict(reads, **{'title': 'Read1 Passed' })
        keys['read1_fail'] = dict(reads, **{'title': 'Read1 Failed' })
        keys['read2_pass'] = dict(reads, **{'title': 'Read2 Passed' })
        keys['read2_fail'] = dict(reads, **{'title': 'Read2 Failed' })
        keys['properly paired_pass'] = dict(reads, **{'title': 'Properly Paired Passed' })
        keys['properly paired_fail'] = dict(reads, **{'title': 'Properly Paired Failed' })
        keys['with itself and mate mapped_pass'] = dict(reads, **{'title': 'Self and mate Passed' })
        keys['with itself and mate mapped_fail'] = dict(reads, **{'title': 'Self and mate Failed' })
        keys['singletons_pass'] = dict(reads, **{'title': 'Mapped Passed' })
        keys['singletons_fail'] = dict(reads, **{'title': 'Mapped Failed' })
        keys['with mate mapped to a different chr_pass'] = dict(reads, **{'title': 'Mate mapped to different chr Passed' })
        keys['with mate mapped to a different chr_fail'] = dict(reads, **{'title': 'Mate mapped to a different chr Failed' })
        keys['with mate mapped to a different chr (mapQ >= 5)_pass'] = dict(reads, **{'title': 'Mate mapped to different chr (mapQ >= 5) Passed' })
        keys['with mate mapped to a different chr (mapQ >= 5)_fail'] = dict(reads, **{'title': 'Mate mapped to a different chr (mapQ >= 5) Failed' })
        
        self.sections.append({
            'name': 'Samtools flagstat output',
            'anchor': 'samtools-flagstat',
            'content': '<p>This module parses the output from <code>samtools flagstat</code>. All numbers in millions.</p>' +
                        plots.beeswarm.plot(self.samtools_flagstat, keys, {'id': 'samtools-flagstat-dp'})
        })
    
    # Return the number of logs that were found
    return len(self.samtools_flagstat)


