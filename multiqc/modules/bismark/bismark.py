#!/usr/bin/env python

""" MultiQC module to parse output from Bismark """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import beeswarm, linegraph, bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

# Log parsing regexes
regexes = {
    'alignment': {
        'total_reads': r"Sequence(?:s| pairs) analysed in total:\s+(\d+)",
        'aligned_reads': r"Number of(?: paired-end)? alignments with a unique best hit(?: from the different alignments)?:\s+(\d+)",
        'no_alignments': r"Sequence(?:s| pairs) with no alignments under any condition:\s+(\d+)",
        'ambig_reads': r"Sequence(?:s| pairs) did not map uniquely:\s+(\d+)",
        'discarded_reads': r"Sequence(?:s| pairs) which were discarded because genomic sequence could not be extracted:\s+(\d+)",
        'total_c': r"Total number of C's analysed:\s+(\d+)",
        'meth_cpg': r"Total methylated C's in CpG context:\s+(\d+)",
        'meth_chg': r"Total methylated C's in CHG context:\s+(\d+)",
        'meth_chh': r"Total methylated C's in CHH context:\s+(\d+)",
        'unmeth_cpg': r"Total unmethylated C's in CpG context:\s+(\d+)",
        'unmeth_chg': r"Total unmethylated C's in CHG context:\s+(\d+)",
        'unmeth_chh': r"Total unmethylated C's in CHH context:\s+(\d+)",
        'percent_cpg_meth': r"C methylated in CpG context:\s+([\d\.]+)%",
        'percent_chg_meth': r"C methylated in CHG context:\s+([\d\.]+)%",
        'percent_chh_meth': r"C methylated in CHH context:\s+([\d\.]+)%",
        'strand_ot': r"CT(?:\/GA)?\/CT:\s+(\d+)\s+\(\(converted\) top strand\)",
        'strand_ctot': r"GA(?:\/CT)?\/CT:\s+(\d+)\s+\(complementary to \(converted\) top strand\)",
        'strand_ctob': r"GA(?:\/CT)?\/GA:\s+(\d+)\s+\(complementary to \(converted\) bottom strand\)",
        'strand_ob': r"CT(?:\/GA)?\/GA:\s+(\d+)\s+\(\(converted\) bottom strand\)",
        'strand_directional': r"Option '--(directional)' specified \(default mode\): alignments to complementary strands \(CTOT, CTOB\) were ignored \(i.e. not performed\)"
    },
    'dedup': {
        'aligned_reads': r"Total number of alignments analysed in .+:\s+(\d+)",
        'dup_reads': r"Total number duplicated alignments removed:\s+(\d+)",
        'dup_reads_percent': r"Total number duplicated alignments removed:\s+\d+\s+\(([\d\.]+)%\)",
        'dedup_reads': r"Total count of deduplicated leftover sequences:\s+(\d+)",
        'dedup_reads_percent': r"Total count of deduplicated leftover sequences:\s+\d+\s+\(([\d\.]+)% of total\)"
    },
    'methextract': {
        'total_c': r"Total number of C's analysed:\s+(\d+)",
        'meth_cpg': r"Total methylated C's in CpG context:\s+(\d+)",
        'meth_chg': r"Total methylated C's in CHG context:\s+(\d+)",
        'meth_chh': r"Total methylated C's in CHH context:\s+(\d+)",
        'unmeth_cpg': r"Total C to T conversions in CpG context:\s+(\d+)",
        'unmeth_chg': r"Total C to T conversions in CHG context:\s+(\d+)",
        'unmeth_chh': r"Total C to T conversions in CHH context:\s+(\d+)",
        'percent_cpg_meth': r"C methylated in CpG context:\s+([\d\.]+)%",
        'percent_chg_meth': r"C methylated in CHG context:\s+([\d\.]+)%",
        'percent_chh_meth': r"C methylated in CHH context:\s+([\d\.]+)%"
    }
}

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Bismark', anchor='bismark',
        href="http://www.bioinformatics.babraham.ac.uk/projects/bismark/",
        info="is a tool to map bisulfite converted sequence reads and determine"\
        " cytosine methylation states.")

        # Set up data structures
        self.bismark_data = {
            'alignment': {},
            'dedup': {},
            'methextract': {},
            'bam2nuc': {}
        }
        self.bismark_mbias_data = {
            'meth': {'CpG_R1' : {}, 'CHG_R1' : {}, 'CHH_R1' : {}, 'CpG_R2' : {}, 'CHG_R2' : {}, 'CHH_R2' : {}},
            'cov': {'CpG_R1' : {}, 'CHG_R1' : {}, 'CHH_R1' : {}, 'CpG_R2' : {}, 'CHG_R2' : {}, 'CHH_R2' : {}}
        }

        # Find and parse bismark alignment reports
        for f in self.find_log_files('bismark/align'):
            parsed_data = self.parse_bismark_report(f['f'], regexes['alignment'])
            if parsed_data is not None:
                # Calculate percent_aligned - doubles as a good check that stuff has worked
                try:
                    parsed_data['percent_aligned'] = (parsed_data['aligned_reads'] / parsed_data['total_reads']) * 100
                except (KeyError, ZeroDivisionError):
                    log.warning('Error calculating percentage for {} - ignoring sample.'.format(f['fn']))
                else:
                    if f['s_name'] in self.bismark_data['alignment']:
                        log.debug("Duplicate alignment sample log found! Overwriting: {}".format(f['s_name']))
                    self.add_data_source(f, section='alignment')
                    self.bismark_data['alignment'][f['s_name']] = parsed_data

        # Find and parse bismark deduplication reports
        for f in self.find_log_files('bismark/dedup'):
            parsed_data = self.parse_bismark_report(f['f'], regexes['dedup'])
            if parsed_data is not None:
                if f['s_name'] in self.bismark_data['dedup']:
                    log.debug("Duplicate deduplication sample log found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f, section='deduplication')
                self.bismark_data['dedup'][f['s_name']] = parsed_data

        # Find and parse bismark methylation extractor reports
        for f in self.find_log_files('bismark/meth_extract'):
            parsed_data = self.parse_bismark_report(f['f'], regexes['methextract'])
            s_name = f['s_name']
            if parsed_data is not None:
                if s_name in self.bismark_data['methextract']:
                    log.debug("Duplicate methylation extraction sample log found! Overwriting: {}".format(s_name))
                self.add_data_source(f, s_name, section='methylation_extraction')
                self.bismark_data['methextract'][s_name] = parsed_data

        # Find and parse M-bias plot data
        for f in self.find_log_files('bismark/m_bias', filehandles=True):
            self.parse_bismark_mbias(f)
            self.add_data_source(f, section='m_bias')

        # Find and parse bam2nuc reports
        for f in self.find_log_files('bismark/bam2nuc', filehandles=True):
            self.parse_bismark_bam2nuc(f)
            self.add_data_source(f, section='bam2nuc')

        # Filters to strip out ignored sample names
        for k in self.bismark_data:
            self.bismark_data[k] = self.ignore_samples(self.bismark_data[k])
        for k in self.bismark_mbias_data['meth']:
            self.bismark_mbias_data['meth'][k] = self.ignore_samples(self.bismark_mbias_data['meth'][k])
        for k in self.bismark_mbias_data['cov']:
            self.bismark_mbias_data['cov'][k] = self.ignore_samples(self.bismark_mbias_data['cov'][k])

        num_parsed = len(self.bismark_data['alignment'])
        num_parsed += len(self.bismark_data['dedup'])
        num_parsed += len(self.bismark_data['methextract'])
        num_parsed += len(self.bismark_mbias_data['meth']['CpG_R1'])
        num_parsed += len(self.bismark_data['bam2nuc'])
        if num_parsed == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Basic Stats Table
        self.bismark_stats_table()

        # Write out to the report
        if len(self.bismark_data['alignment']) > 0:
            self.write_data_file(self.bismark_data['alignment'], 'multiqc_bismark_alignment', sort_cols=True)
            log.info("Found {} bismark alignment reports".format(len(self.bismark_data['alignment'])))
            self.bismark_alignment_chart()

        if len(self.bismark_data['dedup']) > 0:
            self.write_data_file(self.bismark_data['dedup'], 'multiqc_bismark_dedup', sort_cols=True)
            log.info("Found {} bismark dedup reports".format(len(self.bismark_data['dedup'])))
            self.bismark_dedup_chart()

        if len(self.bismark_data['alignment']) > 0:
            self.bismark_strand_chart()

        if len(self.bismark_data['methextract']) > 0:
            self.write_data_file(self.bismark_data['methextract'], 'multiqc_bismark_methextract', sort_cols=True)
            log.info("Found {} bismark methextract reports".format(len(self.bismark_data['methextract'])))
            self.bismark_methlyation_chart()

        if len(self.bismark_mbias_data['meth']['CpG_R1']) > 0:
            self.bismark_mbias_plot()

        if len(self.bismark_data['bam2nuc']) > 0:
            self.write_data_file(self.bismark_data['bam2nuc'], 'multiqc_bismark_bam2nuc', sort_cols=True)
            log.info("Found {} bismark bam2nuc reports".format(len(self.bismark_data['bam2nuc'])))

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

    def parse_bismark_mbias(self, f):
        """ Parse the Bismark M-Bias plot data """
        s = f['s_name']
        self.bismark_mbias_data['meth']['CpG_R1'][s] = {}
        self.bismark_mbias_data['meth']['CHG_R1'][s] = {}
        self.bismark_mbias_data['meth']['CHH_R1'][s] = {}
        self.bismark_mbias_data['cov']['CpG_R1'][s] = {}
        self.bismark_mbias_data['cov']['CHG_R1'][s] = {}
        self.bismark_mbias_data['cov']['CHH_R1'][s] = {}
        self.bismark_mbias_data['meth']['CpG_R2'][s] = {}
        self.bismark_mbias_data['meth']['CHG_R2'][s] = {}
        self.bismark_mbias_data['meth']['CHH_R2'][s] = {}
        self.bismark_mbias_data['cov']['CpG_R2'][s] = {}
        self.bismark_mbias_data['cov']['CHG_R2'][s] = {}
        self.bismark_mbias_data['cov']['CHH_R2'][s] = {}
        key = None
        for l in f['f']:
            if 'context' in l:
                if 'CpG' in l:
                    key = 'CpG'
                elif 'CHG' in l:
                    key = 'CHG'
                elif 'CHH' in l:
                    key = 'CHH'
                if '(R1)' in l:
                    key += '_R1'
                elif '(R2)' in l:
                    key += '_R2'
                else:
                    key += '_R1'
            if key is not None:
                sections = l.split()
                try:
                    pos = int(sections[0])
                    self.bismark_mbias_data['meth'][key][s][pos] = float(sections[3])
                    self.bismark_mbias_data['cov'][key][s][pos] = int(sections[4])
                except (IndexError, ValueError):
                    continue

    def parse_bismark_bam2nuc(self, f):
        """ Parse reports generated by Bismark bam2nuc """
        if f['s_name'] in self.bismark_data['bam2nuc']:
            log.debug("Duplicate deduplication sample log found! Overwriting: {}".format(f['s_name']))
        self.add_data_source(f, section='bam2nuc')
        self.bismark_data['bam2nuc'][f['s_name']] = dict()

        headers = None
        for l in f['f']:
            sections = l.rstrip().split("\t")
            if headers is None:
                headers = sections
            else:
                for i, h in enumerate(headers):
                    if i == 0:
                        k = sections[0]
                    else:
                        key = "{}_{}".format(k, h.lower().replace(' ','_'))
                        self.bismark_data['bam2nuc'][f['s_name']][key] = sections[i]

    def bismark_stats_table(self):
        """ Take the parsed stats from the Bismark reports and add them to the
        basic stats table at the top of the report """

        headers = {
            'alignment': OrderedDict(),
            'dedup': OrderedDict(),
            'methextract': OrderedDict(),
            'bam2nuc': OrderedDict()
        }
        headers['methextract']['percent_cpg_meth'] = {
            'title': '% mCpG',
            'description': '% Cytosines methylated in CpG context',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Greens'
        }
        headers['methextract']['percent_chg_meth'] = {
            'title': '% mCHG',
            'description': '% Cytosines methylated in CHG context',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Oranges'
        }
        headers['methextract']['percent_chh_meth'] = {
            'title': '% mCHH',
            'description': '% Cytosines methylated in CHH context',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Oranges'
        }
        headers['methextract']['total_c'] = {
            'title': "M C's",
            'description': 'Total number of C\'s analysed, in millions',
            'min': 0,
            'scale': 'Purples',
            'modify': lambda x: x / 1000000
        }
        headers['bam2nuc']['C_coverage'] = {
            'title': 'C Coverage',
            'description': 'Cyotosine Coverage',
            'min': 0,
            'suffix': 'X',
            'scale': 'Greens',
            'format': '{:,.2f}'
        }
        headers['dedup']['dup_reads_percent'] = {
            'title': '% Dups',
            'description': 'Percent Duplicated Alignments',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn-rev',
        }
        headers['dedup']['dedup_reads'] = {
            'title': '{} Unique'.format(config.read_count_prefix),
            'description': 'Deduplicated Alignments ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'Greens',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count',
            'hidden': True
        }
        headers['alignment']['aligned_reads'] = {
            'title': '{} Aligned'.format(config.read_count_prefix),
            'description': 'Total Aligned Sequences ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count',
            'hidden': True
        }
        headers['alignment']['percent_aligned'] = {
            'title': '% Aligned',
            'description': 'Percent Aligned Sequences',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }

        self.general_stats_addcols(self.bismark_data['methextract'], headers['methextract'])
        self.general_stats_addcols(self.bismark_data['bam2nuc'], headers['bam2nuc'])
        self.general_stats_addcols(self.bismark_data['dedup'], headers['dedup'])
        self.general_stats_addcols(self.bismark_data['alignment'], headers['alignment'])


    def bismark_alignment_chart (self):
        """ Make the alignment plot """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['aligned_reads']   = { 'color': '#2f7ed8', 'name': 'Aligned Uniquely' }
        keys['ambig_reads']     = { 'color': '#492970', 'name': 'Aligned Ambiguously' }
        keys['no_alignments']   = { 'color': '#0d233a', 'name': 'Did Not Align' }
        keys['discarded_reads'] = { 'color': '#f28f43', 'name': 'No Genomic Sequence' }

        # Config for the plot
        config = {
            'id': 'bismark_alignment',
            'title': 'Bismark Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        self.add_section (
            name = 'Alignment Rates',
            anchor = 'bismark-alignment',
            plot = bargraph.plot(self.bismark_data['alignment'], keys, config)
        )


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
        for sn in self.bismark_data['alignment'].values():
            if 'strand_directional' in sn.keys():
                directional += 1
        if directional == len(self.bismark_data['alignment']):
            keys.pop('strand_ctob', None)
            keys.pop('strand_ctot', None)
            d_mode = 'All samples were run with <code>--directional</code> mode; alignments to complementary strands (CTOT, CTOB) were ignored.'
        elif directional > 0:
            d_mode = '{} samples were run with <code>--directional</code> mode; alignments to complementary strands (CTOT, CTOB) were ignored.'.format(directional)

        # Config for the plot
        config = {
            'id': 'bismark_strand_alignment',
            'title': 'Alignment to Individual Bisulfite Strands',
            'ylab': '% Reads',
            'cpswitch_c_active': False,
            'cpswitch_counts_label': 'Number of Reads'
        }

        self.add_section (
            name = 'Strand Alignment',
            anchor = 'bismark-strands',
            description = d_mode,
            plot = bargraph.plot(self.bismark_data['alignment'], keys, config)
        )


    def bismark_dedup_chart (self):
        """ Make the deduplication plot """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['dedup_reads'] = { 'name': 'Deduplicated reads (remaining)' }
        keys['dup_reads']   = { 'name': 'Duplicate reads (removed)' }

        # Config for the plot
        config = {
            'id': 'bismark_deduplication',
            'title': 'Bismark Deduplication',
            'ylab': '% Reads',
            'cpswitch_c_active': False,
            'cpswitch_counts_label': 'Number of Reads'
        }

        self.add_section (
            name = 'Deduplication',
            anchor = 'bismark-deduplication',
            plot = bargraph.plot(self.bismark_data['dedup'], keys, config)
        )



    def bismark_methlyation_chart (self):
        """ Make the methylation plot """

        # Config for the plot
        keys = OrderedDict()
        defaults = {
            'max': 100,
            'min': 0,
            'suffix': '%',
            'decimalPlaces': 1
        }
        keys['percent_cpg_meth'] = dict(defaults, **{ 'title': 'Methylated CpG' })
        keys['percent_chg_meth'] = dict(defaults, **{ 'title': 'Methylated CHG' })
        keys['percent_chh_meth'] = dict(defaults, **{ 'title': 'Methylated CHH' })

        self.add_section (
            name = 'Cytosine Methylation',
            anchor = 'bismark-methylation',
            plot = beeswarm.plot(self.bismark_data['methextract'], keys, {'id': 'bismark-methylation-dp'})
        )


    def bismark_mbias_plot (self):
        """ Make the M-Bias plot """

        description = '<p>This plot shows the average percentage methylation and coverage across reads. See the \n\
        <a href="https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html#m-bias-plot" target="_blank">bismark user guide</a> \n\
        for more information on how these numbers are generated.</p>'

        pconfig = {
            'id': 'bismark_mbias',
            'title': 'M-Bias',
            'ylab': '% Methylation',
            'xlab': 'Position (bp)',
            'xDecimals': False,
            'ymax': 100,
            'ymin': 0,
            'tt_label': '<b>{point.x} bp</b>: {point.y:.1f}%',
            'data_labels': [
                {'name': 'CpG R1', 'ylab': '% Methylation', 'ymax': 100},
                {'name': 'CHG R1', 'ylab': '% Methylation', 'ymax': 100},
                {'name': 'CHH R1', 'ylab': '% Methylation', 'ymax': 100}
            ]
        }
        datasets = [
            self.bismark_mbias_data['meth']['CpG_R1'],
            self.bismark_mbias_data['meth']['CHG_R1'],
            self.bismark_mbias_data['meth']['CHH_R1']
        ]

        if len(self.bismark_mbias_data['meth']['CpG_R2']) > 0:
            pconfig['data_labels'].append({'name': 'CpG R2', 'ylab': '% Methylation', 'ymax': 100})
            pconfig['data_labels'].append({'name': 'CHG R2', 'ylab': '% Methylation', 'ymax': 100})
            pconfig['data_labels'].append({'name': 'CHH R2', 'ylab': '% Methylation', 'ymax': 100})
            datasets.append(self.bismark_mbias_data['meth']['CpG_R2'])
            datasets.append(self.bismark_mbias_data['meth']['CHG_R2'])
            datasets.append(self.bismark_mbias_data['meth']['CHH_R2'])

        self.add_section (
            name = 'M-Bias',
            anchor = 'bismark-mbias',
            description = description,
            plot = linegraph.plot(datasets, pconfig)
        )
