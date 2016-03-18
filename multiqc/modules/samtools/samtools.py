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
            s_name = f['s_name']
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

        # Add the submodule results to the report output
        self.sections = list()
        self.report_sections()


    def report_sections(self):
        """Create barplots for different stats values"""
        # Genomic Origin Bar Graph
        bases_cats = []
        pairs_cats = []
        reads_cats = ['non-primary alignments']
        if len(self.samtools_data) > 0:
            for s_name in self.samtools_data:
                for metric in self.samtools_data[s_name]:
                    if metric.startswith("reads"):
                        reads_cats.append(metric)
                    elif metric.startswith("bases"):
                        bases_cats.append(metric)
                    elif metric.find("pairs") > -1:
                        pairs_cats.append(metric)
                break
            reads_pconfig = {
                'title': 'Reads statistics',
                'cpswitch': False,
                'cpswitch_c_active': True,
                'stacking': 'normal',
            }
            bases_pconfig= {
                'title': 'Bases statistics',
                'cpswitch': False,
                'cpswitch_c_active': True,
                'stacking': 'normal',
            }
            pairs_pconfig = {
                'title': 'Pairs statistics',
                'cpswitch': False,
                'cpswitch_c_active': True,
                'stacking': 'normal',
            }
            self.sections.append({
                'name': 'Reads mapping statistics',
                'anchor': 'samtools-mapping-statistics',
                'content': self.plot_bargraph(self.samtools_data, reads_cats, reads_pconfig)
            })
            self.sections.append({
                'name': 'Bases mapping statistics',
                'anchor': 'samtools-bases-statistics',
                'content': self.plot_bargraph(self.samtools_data, bases_cats, bases_pconfig)
            })
            # only if data
            if any([any([self.samtools_data[k][v] > 0 for v in pairs_cats]) for k in self.samtools_data]):
                self.sections.append({
                    'name': 'Reads pairs statistics',
                    'anchor': 'samtools-pairs-statistics',
                    'content': self.plot_bargraph(self.samtools_data, pairs_cats, pairs_pconfig)
                })


    def parse_samtools_report (self, raw_data):
        """ Parse the samtools log file. """
        parsed_data = {}
        for line in raw_data.split("\n"):
            if not line.startswith("SN"):
                continue
            fields = line.split("\t")
            parsed_data[fields[1].strip()[:-1]] = float(fields[2].strip())
        if len(parsed_data) == 0: return None
        return parsed_data


    def samtools_stats_table(self):
        """ Take the parsed stats from the samtools report and add them to the
        basic stats table at the top of the report """
        headers = OrderedDict()
        headers['raw total sequences'] = {
            'title': 'M Total seqs',
            'description': 'Total sequences in the bam file',
            'min': 0,
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['reads mapped'] = {
            'title': 'M Reads Mapped',
            'description': 'Reads Mapped in the bam file',
            'min': 0,
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['non-primary alignments'] = {
            'title': 'M Non-Primary Alignments',
            'description': 'Non primary alignment (millions)',
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['error rate'] = {
            'title': 'Error rate',
            'description': 'Error rate using CIGAR',
            'min': 0,
            'max': 100,
            'scale': 'OrRd',
            'format': '{:.2f}%',
            'modify': lambda x: x * 100.0
        }
        self.general_stats_addcols(self.samtools_data, headers)


