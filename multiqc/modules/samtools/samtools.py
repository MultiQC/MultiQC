#!/usr/bin/env python

""" MultiQC module to parse output from samtools stats command"""

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='samtools',
        anchor='samtools', target='samtools',
        href='http://www.htslib.org',
        info="Samtools is a suite of programs for interacting"\
        " with high-throughput sequencing data.")
        self.samtools_data = dict()
        for f in self.find_log_files(config.sp['samtools']):
            s_name = f['s_name'].split(".txt")[0]
            parsed_data = self.parse_samtools_report(f['f'])
            if parsed_data is not None:
                if s_name in self.samtools_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f)
                self.samtools_data[s_name] = parsed_data

        if len(self.samtools_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.samtools_data)))

        # Write parsed report data to a file
        self.write_data_file(self.samtools_data, 'multiqc_samtools')

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.samtools_stats_table()


    def parse_samtools_report (self, raw_data):
        """ Parse the samtools log file. """

        regexes = {
                'raw_total_seqs': r"raw\stotal\ssequences:\s+(\d+)",
            'reads_mapped': r"reads\smapped:\s+(\d+)",
            'non_primary': r"non-primary\salignments:\s+(\d+)",
            'reads_paired': r"reads\smapped\sand\spaired:\s+(\d+)",
            'reads_duplicated': r"reads\sduplicated:\s+(\d+)",
            'error_rate': r"error\srate:\s+(\d+\.\d+e-\d+)\s+#",
            'pairs_on_diff_chrom': r"pairs\son\sdifferent\schromosomes:\s+(\d+)"
        }
        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                value = float(r_search.group(1))
                if value > 0:
                    parsed_data[k] = value
        if len(parsed_data) == 0: return None
        return parsed_data


    def samtools_stats_table(self):
        """ Take the parsed stats from the samtools report and add them to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['raw_total_seqs'] = {
            'title': 'Total seqs',
            'description': 'Total sequences in the bam file',
            'min': 0,
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['raw_total_seqs'] = {
            'title': 'Total seqs',
            'description': 'Total sequences in the bam file',
            'min': 0,
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['error_rate'] = {
            'title': 'Error rate',
            'description': 'Error rate using CIGAR',
            'min': 0,
            'max': 1,
            'scale': 'OrRd',
            'format': '{:.4f}%',
            'modify': lambda x: x * 100.0
        }
        headers['non_primary'] = {
            'title': 'Non primary alignment',
            'description': 'Non primary alignment (millions)',
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
            }
        headers['reads_duplicated'] = {
            'title': 'Duplicatd reads',
            'description': 'PCR or optical duplicate bit set (millions)',
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['reads_paired'] = {
            'title': 'Paired mapped reads',
            'description': 'Paired mapped reads (millions)',
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
            }
        headers['pairs_on_diff_chrom'] = {
            'title': 'Non proper paired *',
            'description': 'Paired in different chromosomes (millions)',
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
            }
        self.general_stats_addcols(self.samtools_data, headers)


