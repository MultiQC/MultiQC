#!/usr/bin/env python

""" MultiQC module to parse output from Bismark """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import io
import json
import logging
import os
import re

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

# Log parsing regexes
regexes = {
    'alignment': {
        'total_reads': r"^Sequence(?:s| pairs) analysed in total:\s+(\d+)$",
        'aligned_reads': r"^Number of(?: paired-end)? alignments with a unique best hit(?: from the different alignments)?:\s+(\d+)$",
        'no_alignments': r"^Sequence(?:s| pairs) with no alignments under any condition:\s+(\d+)$",
        'ambig_reads': r"^Sequence(?:s| pairs) did not map uniquely:\s+(\d+)$",
        'discarded_reads': r"^Sequence(?:s| pairs) which were discarded because genomic sequence could not be extracted:\s+(\d+)$",
        'total_c': r"^Total number of C's analysed:\s+(\d+)$",
        'meth_cpg': r"^Total methylated C's in CpG context:\s+(\d+)",
        'meth_chg': r"^Total methylated C's in CHG context:\s+(\d+)",
        'meth_chh': r"^Total methylated C's in CHH context:\s+(\d+)",
        'unmeth_cpg': r"^Total unmethylated C's in CpG context:\s+(\d+)",
        'unmeth_chg': r"^Total unmethylated C's in CHG context:\s+(\d+)",
        'unmeth_chh': r"^Total unmethylated C's in CHH context:\s+(\d+)",
        'percent_cpg_meth': r"^C methylated in CpG context:\s+([\d\.]+)%",
        'percent_chg_meth': r"^C methylated in CHG context:\s+([\d\.]+)%",
        'percent_chh_meth': r"^C methylated in CHH context:\s+([\d\.]+)%",
        'strand_ot': r"^CT(?:\/GA)?\/CT:\s+(\d+)\s+\(\(converted\) top strand\)$",
        'strand_ctot': r"^GA(?:\/CT)?\/CT:\s+(\d+)\s+\(complementary to \(converted\) top strand\)$",
        'strand_ctob': r"^GA(?:\/CT)?\/GA:\s+(\d+)\s+\(complementary to \(converted\) bottom strand\)$",
        'strand_ob': r"^CT(?:\/GA)?\/GA:\s+(\d+)\s+\(\(converted\) bottom strand\)$",
        'strand_directional': r"^Option '--(directional)' specified \(default mode\): alignments to complementary strands \(CTOT, CTOB\) were ignored \(i.e. not performed\)$"
    },
    'dedup': {
        # 'aligned_reads' overwrites previous, but I trust this more
        'aligned_reads': r"^Total number of alignments analysed in .+:\s+(\d+)$",
        'dup_reads': r"^Total number duplicated alignments removed:\s+(\d+)",
        'dup_reads_percent': r"^Total number duplicated alignments removed:\s+\d+\s+\(([\d\.]+)%\)",
        'dedup_reads': r"^Total count of deduplicated leftover sequences:\s+(\d+)",
        'dedup_reads_percent': r"^Total count of deduplicated leftover sequences:\s+\d+\s+\(([\d\.]+)% of total\)"
    },
    'methextract': {
        # These calls are typically done after deduplication.
        # Overwrites what was found in the alignment report, but now
        # after deduplication (if that was done)
        'total_c': r"^Total number of C's analysed:\s+(\d+)$",
        'meth_cpg': r"^Total methylated C's in CpG context:\s+(\d+)",
        'meth_chg': r"^Total methylated C's in CHG context:\s+(\d+)",
        'meth_chh': r"^Total methylated C's in CHH context:\s+(\d+)",
        'unmeth_cpg': r"^Total C to T conversions in CpG context:\s+(\d+)",
        'unmeth_chg': r"^Total C to T conversions in CHG context:\s+(\d+)",
        'unmeth_chh': r"^Total C to T conversions in CHH context:\s+(\d+)",
        'percent_cpg_meth': r"^C methylated in CpG context:\s+([\d\.]+)%",
        'percent_chg_meth': r"^C methylated in CHG context:\s+([\d\.]+)%",
        'percent_chh_meth': r"^C methylated in CHH context:\s+([\d\.]+)%"
    }
}

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Bismark', anchor='bismark', 
        href="http://www.bioinformatics.babraham.ac.uk/projects/bismark/",
        info="is a tool to map bisulfite converted sequence reads and determine"\
        " cytosine methylation states.")

        # Find and load any Bismark reports
        self.bismark_data = {
            'alignment': {},
            'dedup': {},
            'methextract': {},
            'merged': defaultdict(lambda: dict())
        }
        # Find and parse bismark alignment reports
        for f in self.find_log_files(['_PE_report.txt', '_SE_report.txt']):
            parsed_data = self.parse_bismark_report(f['f'], regexes['alignment'])
            if parsed_data is not None:
                if f['s_name'] in self.bismark_data['alignment']:
                    log.debug("Duplicate alignment sample log found! Overwriting: {}".format(f['s_name']))
                self.bismark_data['alignment'][f['s_name']] = parsed_data
        
        # Find and parse bismark deduplication reports
        for f in self.find_log_files('.deduplication_report.txt'):
            parsed_data = self.parse_bismark_report(f['f'], regexes['dedup'])
            if parsed_data is not None:
                if f['s_name'] in self.bismark_data['dedup']:
                    log.debug("Duplicate deduplication sample log found! Overwriting: {}".format(f['s_name']))
                self.bismark_data['dedup'][f['s_name']] = parsed_data
        
        # Find and parse bismark methylation extractor reports
        for f in self.find_log_files('_splitting_report.txt'):
            parsed_data = self.parse_bismark_report(f['f'], regexes['methextract'])
            if parsed_data is not None:
                if f['s_name'] in self.bismark_data['methextract']:
                    log.debug("Duplicate methylation extraction sample log found! Overwriting: {}".format(f['s_name']))
                self.bismark_data['methextract'][f['s_name']] = parsed_data
        
        # Merge the keys into a merged dict (in order)
        for t in ['alignment','dedup','methextract']:
            for sn in self.bismark_data[t].keys():
                for k, v in self.bismark_data[t][sn].items():
                    self.bismark_data['merged'][sn][k] = v
        
        # Calculate percent_aligned
        for sn in self.bismark_data['merged']:
            aln = self.bismark_data['merged'][sn]['aligned_reads']
            tot = self.bismark_data['merged'][sn]['total_reads']
            self.bismark_data['merged'][sn]['percent_aligned'] = (aln / tot) * 100
        
        if len(self.bismark_data['merged']) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.bismark_data['merged'])))

        # Write parsed report data to a file
        self.write_csv_file(self.bismark_data['merged'], 'multiqc_bismark.txt')

        self.sections = list()
        
        # Basic Stats Table
        self.bismark_stats_table()

        # Section 1 - Column chart of alignment stats
        self.sections.append({
            'name': 'Alignment Rates',
            'anchor': 'bismark-alignment',
            'content': self.bismark_alignment_chart()
        })

        # Section 2 - Methylation percentages
        self.sections.append({
            'name': 'Cytosine Methylation',
            'anchor': 'bismark-methylation',
            'content': self.bismark_methlyation_chart()
        })

        # Section 3 - Strand Alignments
        self.sections.append({
            'name': 'Strand Alignment',
            'anchor': 'bismark-strands',
            'content': self.bismark_strand_chart()
        })

    def parse_bismark_report(self, report, regexes):
        """ Search a bismark report with a set of regexes """
        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, report, re.MULTILINE)
            if r_search:
                try:
                    parsed_data[k] = float(r_search.group(1))
                except ValueError:
                    parsed_data[k] = r_search.group(1) # NaN
        if len(parsed_data) == 0: return None
        return parsed_data

    def bismark_stats_table(self):
        """ Take the parsed stats from the Bismark reports and add them to the
        basic stats table at the top of the report """
        
        headers = OrderedDict()
        headers['percent_cpg_meth'] = {
            'title': '% Meth',
            'description': '% Cytosines methylated in CpG context (alignment)',
            'max': 100,
            'min': 0,
            'scale': 'Greens',
            'format': '{:.1f}%'
        }
        headers['total_c'] = {
            'title': "M C's",
            'description': 'Total number of C\'s analysed, in millions (alignment)',
            'min': 0,
            'scale': 'Purples',
            'modify': lambda x: x / 1000000
        }
        headers['dup_reads_percent'] = {
            'title': '% Dups',
            'description': 'Percent Duplicated Alignments',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn-rev',
            'format': '{:.1f}%'
        }
        headers['dedup_reads'] = {
            'title': 'M Unique',
            'description': 'Deduplicated Alignments (millions)',
            'min': 0,
            'scale': 'Greens',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['aligned_reads'] = {
            'title': 'M Aligned',
            'description': 'Total Aligned Sequences (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['percent_aligned'] = {
            'title': '% Aligned',
            'description': 'Percent Aligned Sequences',
            'max': 100,
            'min': 0,
            'scale': 'YlGn',
            'format': '{:.1f}%',
        }

        self.general_stats_addcols(self.bismark_data['merged'], headers)
        

    def bismark_alignment_chart (self):
        """ Make the alignment plot """
        
        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['dup_reads']       = { 'color': '#8bbc21', 'name': 'Deduplicated Unique Alignments' }
        keys['dedup_reads']     = { 'color': '#2f7ed8', 'name': 'Duplicated Unique Alignments' }
        keys['aligned_reads']   = { 'color': '#2f7ed8', 'name': 'Aligned Uniquely' }
        keys['ambig_reads']     = { 'color': '#492970', 'name': 'Aligned Ambiguously' }
        keys['no_alignments']   = { 'color': '#0d233a', 'name': 'Did Not Align' }
        keys['discarded_reads'] = { 'color': '#f28f43', 'name': 'No Genomic Sequence' }
        
        if len(self.bismark_data['dedup']) == 0:
            keys.pop('dup_reads', None)
            keys.pop('dedup_reads', None)
        else:
            keys.pop('aligned_reads', None)
        
        # Config for the plot
        config = {
            'title': 'Bismark Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return self.plot_bargraph(self.bismark_data['merged'], keys, config)


    def bismark_methlyation_chart (self):
        """ Make the methylation plot """
        
        # Specify the order of the different possible categories
        cats = [OrderedDict(), OrderedDict(), OrderedDict()]
        cats[0]['meth_cpg'] =   {'color': '#0d233a', 'name': 'Methylated CpG'}
        cats[0]['unmeth_cpg'] = {'color': '#2f7ed8', 'name': 'Unmethylated CpG'}
        cats[1]['meth_chg'] =   {'color': '#1aadce', 'name': 'Methylated CHG'}
        cats[1]['unmeth_chg'] = {'color': '#8bbc21', 'name': 'Unmethylated CHG'}
        cats[2]['meth_chh'] =   {'color': '#492970', 'name': 'Methylated CHH'}
        cats[2]['unmeth_chh'] = {'color': '#910000', 'name': 'Unmethylated CHH'}
        
        # Config for the plot
        config = {
            'title': 'Cytosine Methylation',
            'ylab': '% Calls',
            'cpswitch_c_active': False,
            'cpswitch_counts_label': 'Number of Calls',
            'data_labels': ['CpG', 'CHG', 'CHH']
        }
        
        # Need to supply three data dicts
        data = [self.bismark_data['merged'], self.bismark_data['merged'], self.bismark_data['merged']]
        
        return self.plot_bargraph(data, cats, config)

    def bismark_strand_chart (self):
        """ Make the strand alignment plot """
        
        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['strand_ob']   = { 'name': 'Original bottom strand' }
        keys['strand_ctob'] = { 'name': 'Complementary to original bottom strand' }
        keys['strand_ctot'] = { 'name': 'Complementary to original top strand' }
        keys['strand_ot']   = { 'name': 'Original top strand' }
        
        # See if we have any directional samples
        directional = 0
        d_mode = ''
        for sn in self.bismark_data['merged'].values():
            if 'strand_directional' in sn.keys():
                directional += 1
        if directional == len(self.bismark_data['merged']):
            keys.pop('strand_ctob', None)
            keys.pop('strand_ctot', None)
            d_mode = '<p>All samples were run with <code>--directional</code> mode; alignments to complementary strands (CTOT, CTOB) were ignored.</p>'
        elif directional > 0:
            d_mode = '<p>{} samples were run with <code>--directional</code> mode; alignments to complementary strands (CTOT, CTOB) were ignored.</p>'.format(directional)
        
        # Config for the plot
        config = {
            'title': 'Alignment to Individual Bisulfite Strands',
            'ylab': '% Reads',
            'cpswitch_c_active': False,
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return d_mode + self.plot_bargraph(self.bismark_data['merged'], keys, config)
