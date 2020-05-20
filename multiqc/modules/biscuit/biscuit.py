#!/usr/bin/env python

''' MultiQC module to parse output from BISCUITqc '''

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc.plots import linegraph, bargraph, table, beeswarm
from multiqc.modules.base_module import BaseMultiqcModule

# Initialize the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    '''
    Inherits from base module class
    Initialize data structures and prepare BISCUIT report
    Inputs:
        No inputs
    Returns:
        BISCUIT report for MultiQC
    '''
    def __init__(self):
        # Initialize the parent object
        super(MultiqcModule, self).__init__(name='BISCUIT', anchor='biscuit',
            href='https://github.com/huishenlab/biscuit',
            info='is a tool to map bisulfite converted DNA sequence reads and' \
            ' determine cytosine methylation states.')

        # Set up data structures
        self.mdata = {
            # General statistics
            'align_mapq': {},
            'align_strand': {},
            'align_isize': {},
            # Duplicate reporting
            'dup_report': {},
            # Uniformity
            'qc_cv': {},
            # Base coverage
            'covdist_all_base_botgc': {},
            'covdist_all_base': {},
            'covdist_all_base_topgc': {},
            'covdist_q40_base_botgc': {},
            'covdist_q40_base': {},
            'covdist_q40_base_topgc': {},
            # CpG coverage
            'covdist_all_cpg_botgc': {},
            'covdist_all_cpg': {},
            'covdist_all_cpg_topgc': {},
            'covdist_q40_cpg_botgc': {},
            'covdist_q40_cpg': {},
            'covdist_q40_cpg_topgc': {},
            # Cytosine retention
            'cpg_retention_readpos': {},
            'cph_retention_readpos': {},
            'base_avg_retention_rate': {},
            'read_avg_retention_rate': {}
        }

        file_suffixes = [
            # General statistics
            '_mapq_table',
            '_strand_table',
            '_isize_table',
            # Duplicate reporting
            '_dup_report',
            # Uniformity
            '_cv_table',
            # Base coverage
            '_covdist_all_base_botgc_table',
            '_covdist_all_base_table',
            '_covdist_all_base_topgc_table',
            '_covdist_q40_base_botgc_table',
            '_covdist_q40_base_table',
            '_covdist_q40_base_topgc_table',
            # CpG coverage
            '_covdist_all_cpg_botgc_table',
            '_covdist_all_cpg_table',
            '_covdist_all_cpg_topgc_table',
            '_covdist_q40_cpg_botgc_table',
            '_covdist_q40_cpg_table',
            '_covdist_q40_cpg_topgc_table',
            # Cytosine retention
            '_CpGRetentionByReadPos',
            '_CpHRetentionByReadPos',
            '_totalBaseConversionRate',
            '_totalReadConversionRate'
        ]

        # Find and parse alignment reports
        for k in self.mdata:
            for f in self.find_log_files('biscuit/{}'.format(k)):
                # Add source file to multiqc_sources.txt
                self.add_data_source(f)

                # Clean s_name before further processing
                s_name = self.clean_s_name(f['s_name'], f['root'])

                # Clean file suffixes unique to biscuit
                for suffix in file_suffixes:
                    s_name = s_name.replace(suffix, '')

                if s_name in self.mdata[k]:
                    log.debug('Duplicate sample name found in {}! Overwriting: {}'.format(f['fn'], s_name))

                self.mdata[k][s_name] = getattr(self, 'parse_logs_{}'.format(k))(f['f'], f['fn'])

        for k in self.mdata:
            self.mdata[k] = self.ignore_samples(self.mdata[k])

        n_reports = sum([len(self.mdata[k]) for k in self.mdata])

        if n_reports == 0:
            raise UserWarning

        # Basic stats table
        self.biscuit_stats_table()

        # Write out BISCUIT MultiQC report
        log.info('Found {} reports'.format(n_reports))
        for k in self.mdata:
            if len(self.mdata[k]) > 0:
                log.debug('Found {} {} reports'.format(len(self.mdata[k]), k))
                self.write_data_file(self.mdata[k],
                                     'multiqc_biscuit_{}'.format(k),
                                     data_format='json')
                getattr(self, 'chart_{}'.format(k))()

    def biscuit_stats_table(self):
        '''
        Create general statistics table for BISCUIT data
        Inputs:
            Uses mdata['align_mapq'] and mdata['dup_report']
        Returns:
            Add columns to MultiQC general statistics table
        '''
        pd = {}

        for sid, dd in self.mdata['align_mapq'].items():
            if len(dd) > 0:
                total = sum([int(_) for _ in dd.values()])
                pd[sid] = {'aligned': 100.0 * float(total - int(dd['unmapped'])) / total}

        for sid, dd in self.mdata['dup_report'].items():
            if sid not in pd:
                pd[sid] = {}
            if 'all' in dd and dd['all'] != -1:
                pd[sid]['dup_all'] = dd['all']
            if 'q40' in dd and dd['q40'] != -1:
                pd[sid]['dup_q40'] = dd['q40']

        pheader = OrderedDict()
        pheader['aligned'] = {'title': '% Aligned', 'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Reds'}
        pheader['dup_all'] = {'title': 'Dup. % for All Reads', 'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Reds'}
        pheader['dup_q40'] = {'title': 'Dup. % for Q40 Reads', 'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Purples'}

        self.general_stats_addcols(pd, pheader)

    ########################################
    #####  General Mapping Information #####
    ########################################
    def parse_logs_align_mapq(self, f, fn):
        '''
        Parse _mapq_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of aligned mapq data
        '''
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug('No data available in {}. Will not fill corresponding entries.'.format(fn))
            return {}

        data = {}
        for l in file_data:
            s = l.split()
            data[s[0]] = s[1] # data[MAPQ] = number of reads

        return data

    def chart_align_mapq(self):
        '''
        Chart _mapq_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Mapping Overview and Mapping Quality
            Distribution charts
        '''

        # Mapping Overview
        pd = {}

        for sid, dd in self.mdata['align_mapq'].items():
            if len(dd) > 0:
                pd[sid] = {'opt_align': 0, 'sub_align': 0, 'not_align': 0}
                for mapq, cnt in dd.items():
                    if mapq == 'unmapped':
                        pd[sid]['not_align'] += int(cnt)
                    elif int(mapq) >= 40:
                        pd[sid]['opt_align'] += int(cnt)
                    else:
                        pd[sid]['sub_align'] += int(cnt)

        pheader = OrderedDict()
        pheader['opt_align'] = {'color': '#a6cee3', 'name': 'Optimally Aligned Reads'}
        pheader['sub_align'] = {'color': '#1f78b4', 'name': 'Suboptimally Aligned Reads'}
        pheader['not_align'] = {'color': '#b2df8a', 'name': 'Unaligned Reads'}

        pconfig = {
            'id': 'biscuit_mapping_overview',
            'title': 'BISCUIT: Mapping Overview',
            'ylab': 'Number of Reads',
            'cpswitch_counts_label': '# Reads'
        }

        self.add_section(
            name = 'Mapping Overview',
            anchor = 'biscuit-mapping-overview',
            description = 'For primary alignments, shows the number of ' \
            'optimally aligned reads (defined by MAPQ>=40), suboptimally ' \
            'aligned reads (MAPQ<40), and unmapped reads. See Help for more ' \
            'details.',
            helptext = 'A good library should have a high fraction of reads ' \
            'that are optimally aligned. Note, suboptimally aligned reads ' \
            'include both non-unique alignments and imperfect alignments.',
            plot = bargraph.plot(pd, pheader, pconfig)
        )

        # Mapping Quality Distribution
        total = {}
        for sid, dd in self.mdata['align_mapq'].items():
            if len(dd) > 0:
                total[sid] = sum([int(cnt) for _, cnt in dd.items() if _ != 'unmapped'])

        pd_mapq = {}
        for sid, dd in self.mdata['align_mapq'].items():
            if len(dd) > 0:
                cnts = []
                for mapq in range(61):
                    if str(mapq) in dd:
                        cnts.append(100.0 * float(dd[str(mapq)]) / total[sid])
                    else:
                        cnts.append(0)
                pd_mapq[sid] = dict(zip(range(61), cnts))

        pconfig = {
            'id': 'biscuit_mapq',
            'title': 'BISCUIT: Distribution of Mapping Qualities',
            'ymin': 0,
            'xmin': 0,
            'yLabelFormat': '{value}%',
            'tt_label': '<strong>Q{point.x}:</strong> {point.y:.2f}% of mapped reads',
            'ylab': '% of Primary Mapped Reads',
            'xlab': 'Mapping Quality Score'
        }

        self.add_section(
            name = 'Mapping Quality Distribution',
            anchor = 'biscuit-mapq',
            description = 'Shows the percentage of the total number of mapped ' \
            'reads each mapping quality score has (for primary alignments only).',
            plot = linegraph.plot(pd_mapq, pconfig)
        )

    def parse_logs_align_strand(self, f, fn):
        '''
        Parse _strand_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of strand data for reads 1 and 2
        '''
        patterns = [
            r'     R1 \(f\)\:   (\d+)(\s+)(\d+)(\s+)',
            r'     R1 \(r\)\:   (\d+)(\s+)(\d+)(\s+)',
            r'     R2 \(f\)\:   (\d+)(\s+)(\d+)(\s+)',
            r'     R2 \(r\)\:   (\d+)(\s+)(\d+)(\s+)'
        ]

        data = {'read1': {}, 'read2': {}}
        for pat in patterns:
            m = re.search(pat, f, re.MULTILINE)
            if m is not None:
                if (m.group(0)[5:7]) == 'R1':
                    if (m.group(0)[9]) == 'f':
                        data['read1']['ff'] = int(float(m.group(1)))
                        data['read1']['fr'] = int(float(m.group(3)))
                    else:
                        data['read1']['rf'] = int(float(m.group(1)))
                        data['read1']['rr'] = int(float(m.group(3)))
                else:
                    if (m.group(0)[9]) == 'f':
                        data['read2']['ff'] = int(float(m.group(1)))
                        data['read2']['fr'] = int(float(m.group(3)))
                    else:
                        data['read2']['rf'] = int(float(m.group(1)))
                        data['read2']['rr'] = int(float(m.group(3)))

        return data

    def chart_align_strand(self):
        '''
        Chart _strand_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Mapping Strand Distribution chart
        '''

        pd1 = {}
        pd2 = {}
        for sid, dd in self.mdata['align_strand'].items():
            pd1[sid] = dd['read1']
            pd2[sid] = dd['read2']

        pd = [pd1, pd2]

        pheader = [
            OrderedDict([('ff', {'color': '#F53855', 'name': 'ff: Waston-Aligned, Waston-Bisulfite Conversion'}),
                         ('fr', {'color': '#E37B40', 'name': 'fr: Waston-Aligned, Crick-Bisulfite Conversion' }),
                         ('rf', {'color': '#46B29D', 'name': 'rf: Crick-Aligned, Waston-Bisulfite Conversion' }),
                         ('rr', {'color': '#324D5C', 'name': 'rr: Crick-Aligned, Crick-Bisulfite Conversion'  })]),
            OrderedDict([('ff', {'color': '#F53855', 'name': 'ff: Waston-Aligned, Waston-Bisulfite Conversion'}),
                         ('fr', {'color': '#E37B40', 'name': 'fr: Waston-Aligned, Crick-Bisulfite Conversion' }),
                         ('rf', {'color': '#46B29D', 'name': 'rf: Crick-Aligned, Waston-Bisulfite Conversion' }),
                         ('rr', {'color': '#324D5C', 'name': 'rr: Crick-Aligned, Crick-Bisulfite Conversion'  })])
        ]

        pconfig = {
            'id': 'biscuit_strands',
            'title': 'BISCUIT: Mapping Strand Distribution',
            'ylab': 'Number of Reads',
            'cpswitch_counts_label': '# Reads',
            'data_labels': [{'name': 'Read 1'}, {'name': 'Read 2'}]
        }

        # TODO: When PBAT mode is implemented, add comment in help text about
        #       how to interpret PBAT mode results
        self.add_section(
            name = 'Mapping Strand Distribution',
            anchor = 'biscuit-strands',
            description = 'For primary alignments, shows the number of reads ' \
            'mapped to each strand. See Help for more details',
            helptext = 'Most bisulfite libraries map read 1 to the parent ' \
            'strand (`ff` or `rr`) and read 2 to the daughter/synthesized ' \
            'strand (`fr` or `rf`). PBAT and most single-cell/low input ' \
            'libraries often do not follow this assumption.',
            plot = bargraph.plot(pd, pheader, pconfig)
        )

    def parse_logs_align_isize(self, f, fn):
        '''
        Parse _isize_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of insert size data
        '''
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug('No data available in {}. Will not fill corresponding entries.'.format(fn))
            return {'no_data_available': 1}

        data = {'percent': {}, 'readcnt': {}}
        for l in file_data:
            fields = l.split('\t')
            data['percent'][int(fields[0])] = 100.0 * float(fields[1])
            data['readcnt'][int(fields[0])] = float(fields[2]) / 1000000.0

        return data

    def chart_align_isize(self):
        '''
        Chart _isize_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Insert Size Distribution chart
        '''

        pd_p = {}
        pd_r = {}
        for sid, dd in self.mdata['align_isize'].items():
            if 'no_data_available' not in dd.keys():
                pd_p[sid] = dd['percent']
                pd_r[sid] = dd['readcnt']

        pd = [pd_p, pd_r]

        pconfig = {
            'id': 'biscuit_isize',
            'title': 'BISCUIT: Insert Size Distribution',
            'ymin': 0,
            'xmin': 0,
            'yLabelFormat': '{value}',
            'smooth_points': 1000, # limit number of points / smooth data
            'tt_label': '<strong>IS{point.x}:</strong> {point.y:.2f}',
            'xlab': 'Insert Size',
            'ylab': '% of Mapped Reads',
            'data_labels': [{'name': 'Percent of Reads', 'ylab':'% of Mapped Reads'},
                            {'name': 'Millions of Reads', 'ylab':'Millions of Mapped Reads'}]
        }

        self.add_section(
            name = 'Insert Size Distribution',
            anchor = 'biscuit-isize',
            description = 'Shows the distribution of insert sizes. See Help ' \
            'for more details.',
            helptext = 'Insert size is defined as: `(right-most coordinate ' \
            'of reverse-mate read) - (left-most coordinate of forward-mate ' \
            'read)`. Insert sizes are calculated for reads with a "mapped in ' \
            'proper pair" `samtools` flag, and MAPQ >= 40.',
            plot = linegraph.plot(pd, pconfig)
        )

    ########################################
    ####        Duplicate Report        ####
    ########################################
    def parse_logs_dup_report(self, f, fn):
        '''
        Parses _dup_report.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of duplicate fractions
        '''
        patterns = [
            (r'Number of duplicate reads:\t(\d+)',
             r'Number of reads:\t(\d+)',
             'all'),
            (r'Number of duplicate q40-reads:\t(\d+)',
             r'Number of q40-reads:\t(\d+)',
             'q40')
        ]

        data = {}
        for pat_dup, pat_tot, key in patterns:
            m1 = re.search(pat_dup, f, re.MULTILINE)
            m2 = re.search(pat_tot, f, re.MULTILINE)
            if m1 is not None and m2 is not None:
                data[key] = (100.0 * float(m1.group(1)) / float(m2.group(1)))
            else:
                data[key] = -1

        return data

    def chart_dup_report(self):
        '''
        Charts _dup_report.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Duplicate Rates chart
        '''

        pd = dict([(sid, dd) for sid, dd in self.mdata['dup_report'].items() if dd['all'] != -1])

        pheader = OrderedDict()
        pheader['all'] = {'title': 'Overall', 'suffix': '%', 'max': 100, 'min': 0, 'scale': 'Reds'}
        pheader['q40'] = {'title': 'MAPQ >= 40', 'suffix': '%', 'max': 100, 'min': 0, 'scale': 'Purples'}

        pconfig = {
            'id': 'biscuit_dup_report',
            'table_title': 'BISCUIT: Duplicate Rates',
            'save_file': True,
            'sortRows': False
        }

        self.add_section(
            name = 'Duplicate Rates',
            anchor = 'biscuit-dup-report',
            description = 'Shows the percentage of reads that are duplicates ' \
            'out of the total number of reads. See Help for more details.',
            helptext = 'MAPQ >= 40 shows the duplicate rate for reads having ' \
            'a MAPQ >= 40.',
            plot = table.plot(pd, pheader, pconfig)
        )

    ########################################
    ####      Depths and Uniformity     ####
    ########################################
    def parse_logs_qc_cv(self, f, fn):
        '''
        Parses _cv_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of depth uniformity measures
        '''

        data = {}
        targets = ['all_base', 'all_cpg',
                   'q40_base', 'q40_cpg',
                   'all_base_botgc', 'all_cpg_botgc',
                   'q40_base_botgc', 'q40_cpg_botgc',
                   'all_base_topgc', 'all_cpg_topgc',
                   'q40_base_topgc', 'q40_cpg_topgc']
        for t in targets:
            m = re.search('{}\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)'.format(t),
                          f, re.MULTILINE)
            if m is not None:
                data[t] = {'mu': float(m.group(1)),
                           'sigma': float(m.group(2)),
                           'cv': float(m.group(3))}
            else:
                data[t] = {'mu': -1, 'sigma': -1, 'cv': -1}

        return data

    def chart_qc_cv(self):
        '''
        Charts _cv_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Sequencing Depth - Whole Genome chart
        '''

        cats = [('all_base', 'a_b'), ('q40_base', 'q_b'),
                ('all_base_botgc', 'a_b_b'), ('q40_base_botgc', 'q_b_b'),
                ('all_base_topgc', 'a_b_t'), ('q40_base_topgc', 'q_b_t'),
                ('all_cpg', 'a_c'), ('q40_cpg', 'q_c'),
                ('all_cpg_botgc', 'a_c_b'), ('q40_cpg_botgc', 'q_c_b'),
                ('all_cpg_topgc', 'a_c_t'), ('q40_cpg_topgc', 'q_c_t')]

        pd = OrderedDict()
        for sid, dd in self.mdata['qc_cv'].items():
            pd[sid] = OrderedDict()
            for cat, key in cats:
                if cat in dd:
                    if dd[cat]['mu'] != -1:
                        pd[sid]['mu_'+key] = dd[cat]['mu']
                        pd[sid]['cv_'+key] = dd[cat]['cv']

        shared_mean = {'min': 0, 'format': '{:,3f}', 'minRange': 10}
        shared_cofv = {'min': 0, 'format': '{:,3f}', 'minRange': 50}

        pheader = OrderedDict()
        pheader['mu_a_b'] = dict(shared_mean, **{'title': 'All Genome Mean', 'description': 'Mean Sequencing Depth for All Reads'})
        pheader['mu_q_b'] = dict(shared_mean, **{'title': 'Q40 Genome Mean', 'description': 'Mean Sequencing Depth for Q40 Reads'})
        pheader['mu_a_b_b'] = dict(shared_mean, **{'title': 'Low GC All Gen. Mean', 'description': 'Mean Sequencing Depth for All Reads in Low GC-Content Regions'})
        pheader['mu_q_b_b'] = dict(shared_mean, **{'title': 'Low GC Q40 Gen. Mean', 'description': 'Mean Sequencing Depth for Q40 Reads in Low GC-Content Regions'})
        pheader['mu_a_b_t'] = dict(shared_mean, **{'title': 'High GC All Gen. Mean', 'description': 'Mean Sequencing Depth for All Reads in High GC-Content Regions'})
        pheader['mu_q_b_t'] = dict(shared_mean, **{'title': 'High GC Q40 Gen. Mean', 'description': 'Mean Sequencing Depth for Q40 Reads in High GC-Content Regions'})
        pheader['cv_a_b'] = dict(shared_cofv, **{'title': 'All Genome CoV', 'description': 'Sequencing Depth CoV for All Reads'})
        pheader['cv_q_b'] = dict(shared_cofv, **{'title': 'Q40 Genome CoV', 'description': 'Sequencing Depth CoV for Q40 Reads'})
        pheader['cv_a_b_b'] = dict(shared_cofv, **{'title': 'Low GC All Gen. CoV', 'description': 'Sequencing Depth CoV for All Reads in Low GC-Content Regions'})
        pheader['cv_q_b_b'] = dict(shared_cofv, **{'title': 'Low GC Q40 Gen. CoV', 'description': 'Sequencing Depth CoV for Q40 Reads in Low GC-Content Regions'})
        pheader['cv_a_b_t'] = dict(shared_cofv, **{'title': 'High GC All Gen. CoV', 'description': 'Sequencing Depth CoV for All Reads in High GC-Content Regions'})
        pheader['cv_q_b_t'] = dict(shared_cofv, **{'title': 'High GC Q40 Gen. CoV', 'description': 'Sequencing Depth CoV for Q40 Reads in High GC-Content Regions'})

        pheader['mu_a_c'] = dict(shared_mean, **{'title': 'All CpGs Mean', 'description': 'Mean Sequencing Depth for All CpGs'})
        pheader['mu_q_c'] = dict(shared_mean, **{'title': 'Q40 CpGs Mean', 'description': 'Mean Sequencing Depth for Q40 CpGs'})
        pheader['mu_a_c_b'] = dict(shared_mean, **{'title': 'Low GC All CpGs Mean', 'description': 'Mean Sequencing Depth for All CpGs in Low GC-Content Regions'})
        pheader['mu_q_c_b'] = dict(shared_mean, **{'title': 'Low GC Q40 CpGs Mean', 'description': 'Mean Sequencing Depth for Q40 CpGs in Low GC-Content Regions'})
        pheader['mu_a_c_t'] = dict(shared_mean, **{'title': 'High GC All CpGs Mean', 'description': 'Mean Sequencing Depth for All CpGs in High GC-Content Regions'})
        pheader['mu_q_c_t'] = dict(shared_mean, **{'title': 'High GC Q40 CpGs Mean', 'description': 'Mean Sequencing Depth for Q40 CpGs in High GC-Content Regions'})
        pheader['cv_a_c'] = dict(shared_cofv, **{'title': 'All CpGs CoV', 'description': 'Sequencing Depth CoV for All CpGs'})
        pheader['cv_q_c'] = dict(shared_cofv, **{'title': 'Q40 CpGs CoV', 'description': 'Sequencing Depth CoV for Q40 CpGs'})
        pheader['cv_a_c_b'] = dict(shared_cofv, **{'title': 'Low GC All CpGs CoV', 'description': 'Sequencing Depth CoV for All CpGs in Low GC-Content Regions'})
        pheader['cv_q_c_b'] = dict(shared_cofv, **{'title': 'Low GC Q40 CpGs CoV', 'description': 'Sequencing Depth CoV for Q40 CpGs in Low GC-Content Regions'})
        pheader['cv_a_c_t'] = dict(shared_cofv, **{'title': 'High GC All CpGs CoV', 'description': 'Sequencing Depth CoV for All CpGs in High GC-Content Regions'})
        pheader['cv_q_c_t'] = dict(shared_cofv, **{'title': 'High GC Q40 CpGs CoV', 'description': 'Sequencing Depth CoV for Q40 CpGs in High GC-Content Regions'})

        pconfig = {
            'id': 'biscuit_seq_depth',
            'table_title': 'BISCUIT: Sequencing Depth',
            'save_file': True,
            'sortRows': False
        }

        self.add_section(
            name = 'Sequencing Depth Statistics',
            anchor = 'biscuit-seq-depth',
            description = 'Shows the sequence depth mean and uniformity ' \
            'measured by the Coefficient of Variation (CoV, defined as ' \
            'std. dev./mean). See Help for more details.',
            helptext = 'The "Genome" (Gen.) show statistics for all bases ' \
            'across the entire genome, while "CpGs" shows the corresponding ' \
            'statistics for CpGs. "All" shows statistics for any mapped ' \
            'bases/CpGs, while "Q40" shows statistics only those bases/CpGs ' \
            'with MAPQ >= 40. "High GC" and "low GC" shows bases/CpGs that ' \
            'overlap with the top and bottom 10% of 100bp windows for ' \
            'GC-content, respectively.',
            plot = beeswarm.plot(pd, pheader, pconfig)
        )

    ########################################
    #### Base Coverage and CpG Coverage ####
    ########################################
    def parse_logs_covdist_all_base(self, f, fn):
        '''
        Parses _covdist_all_base_botgc_table.txt
               _covdist_all_base_table.txt
               _covdist_all_base_topgc_table.txt
               _covdist_all_cpg_botgc_table.txt
               _covdist_all_cpg_table.txt
               _covdist_all_cpg_topgc_table.txt
               _covdist_q40_base_botgc_table.txt
               _covdist_q40_base_table.txt
               _covdist_q40_base_topgc_table.txt
               _covdist_q40_cpg_botgc_table.txt
               _covdist_q40_cpg_table.txt
               _covdist_q40_cpg_topgc_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of coverage distributions up to 30X data
        '''
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug('No data available in {}. Will not fill corresponding entries.'.format(fn))
            return dict(zip([i for i in range(31)], [-1 for _ in range(31)]))

        dd = {}
        for l in file_data:
            fields = l.split()
            dd[int(float(fields[0]))] = int(float(fields[1]))

        covs = sorted([k for k in dd])[:31]
        _ccov_cnt = sum(dd.values())

        ccov_cnts = []
        for cov in covs:
            ccov_cnts.append(_ccov_cnt/1000000.0)
            _ccov_cnt -= dd[cov]

        return dict(zip(covs, ccov_cnts))

    def parse_logs_covdist_all_base_botgc(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_all_base_topgc(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_base(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_base_botgc(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_base_topgc(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_all_cpg(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_all_cpg_botgc(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_all_cpg_topgc(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_cpg(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_cpg_botgc(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_cpg_topgc(self, f, fn):
        ''' Handled by parse_logs_covdist_all_base() '''
        return self.parse_logs_covdist_all_base(f, fn)

    def chart_covdist_all_base(self):
        '''
        Charts _covdist_all_base_botgc_table.txt
               _covdist_all_base_table.txt
               _covdist_all_base_topgc_table.txt
               _covdist_all_cpg_botgc_table.txt
               _covdist_all_cpg_table.txt
               _covdist_all_cpg_topgc_table.txt
               _covdist_q40_base_botgc_table.txt
               _covdist_q40_base_table.txt
               _covdist_q40_base_topgc_table.txt
               _covdist_q40_cpg_botgc_table.txt
               _covdist_q40_cpg_table.txt
               _covdist_q40_cpg_topgc_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Cumulative Coverage chart
        '''

        pd = [
            self.mdata['covdist_all_base'],
            self.mdata['covdist_q40_base'],
            self.mdata['covdist_all_cpg'],
            self.mdata['covdist_q40_cpg'],
            self.mdata['covdist_all_base_botgc'],
            self.mdata['covdist_q40_base_botgc'],
            self.mdata['covdist_all_cpg_botgc'],
            self.mdata['covdist_q40_cpg_botgc'],
            self.mdata['covdist_all_base_topgc'],
            self.mdata['covdist_q40_base_topgc'],
            self.mdata['covdist_all_cpg_topgc'],
            self.mdata['covdist_q40_cpg_topgc']
        ]

        pconfig = {
            'id': 'biscuit_cumulative',
            'title': 'BISCUIT: Cumulative Coverage',
            'ymin': 0,
            'xLabelFormat': '{value}X',
            'tt_label': '<strong>{point.x}X:</strong> {point.y:.2f}M',
            'xlab': 'Coverage',
            'ylab': 'Millions of Bases',
            'data_labels': [{'name': 'All Bases', 'ylab': 'Millions of Bases'},
                            {'name': 'Q40 Bases', 'ylab': 'Millions of Bases'},
                            {'name': 'All CpGs', 'ylab': 'Millions of CpGs'},
                            {'name': 'Q40 CpGs', 'ylab': 'Millions of CpGs'},
                            {'name': 'Low GC All Bases', 'ylab': 'Millions of Bases'},
                            {'name': 'Low GC Q40 Bases', 'ylab': 'Millions of Bases'},
                            {'name': 'Low GC All CpGs', 'ylab': 'Millions of CpGs'},
                            {'name': 'Low GC Q40 CpGs', 'ylab': 'Millions of CpGs'},
                            {'name': 'High GC All Bases', 'ylab': 'Millions of Bases'},
                            {'name': 'High GC Q40 Bases', 'ylab': 'Millions of Bases'},
                            {'name': 'High GC All CpGs', 'ylab': 'Millions of CpGs'},
                            {'name': 'High GC Q40 CpGs', 'ylab': 'Millions of CpGs'}]
        }

        self.add_section(
            name = 'Cumulative Coverage',
            anchor = 'biscuit-cumulative-coverage',
            description = 'Shows the number of bases or CpGs covered by a given ' \
            'number of reads. See Help for more details.',
            helptext = '"All" shows the coverage for any mapped reads, while ' \
            '"Q40" shows only those reads with MAPQ >= 40. "High GC" and ' \
            '"low GC" shows reads that overlap with the top and bottom 10% ' \
            'of 100bp windows for GC-content, respectively.',
            plot = linegraph.plot(pd, pconfig)
        )

    def chart_covdist_all_base_botgc(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_all_base_topgc(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_q40_base(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_q40_base_botgc(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_q40_base_topgc(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_all_cpg(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_all_cpg_botgc(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_all_cpg_topgc(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_q40_cpg(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_q40_cpg_botgc(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    def chart_covdist_q40_cpg_topgc(self):
        ''' Handled by chart_covdist_all_base() '''
        pass

    ########################################
    ####          CpG Retention         ####
    ########################################
    def parse_logs_cpg_retention_readpos(self, f, fn):
        '''
        Parses _CpGRetentionByReadPos.txt
               _CpHRetentionByReadPos.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of fraction of retained cytosines for reads 1 and 2
            in either a CpH or CpG context
        '''
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug('No data available in {}. Will not fill corresponding entries.'.format(fn))
            return {'no_data_available': 1}

        r1 = {'C': {}, 'R': {}}
        r2 = {'C': {}, 'R': {}}
        for l in file_data:
            fields = l.strip().split('\t')

            if fields[0] not in ['1', '2'] or fields[2] not in ['C', 'R']:
                return {}
            if fields[0] == '1':
                r1[fields[2]][int(fields[1])] = int(fields[3])
            elif fields[0] == '2':
                r2[fields[2]][int(fields[1])] = int(fields[3])

        r1rate = OrderedDict()
        for k in sorted(r1['C'].keys()):
            if k in r1['R']:
                r1rate[k] = 100.0 * float(r1['R'][k]) / (r1['R'][k] + r1['C'][k])

        r2rate = OrderedDict()
        for k in sorted(r2['C'].keys()):
            if k in r2['R']:
                r2rate[k] = 100.0 * float(r2['R'][k]) / (r2['R'][k] + r2['C'][k])

        return {'1': r1rate, '2': r2rate}

    def chart_cpg_retention_readpos(self):
        '''
        Charts _CpGRetentionByReadPos.txt
               _CpHRetentionByReadPos.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Retenion vs. Base Position in Read chart
        '''

        pd = [dict([(sid, dd['1']) for sid, dd in self.mdata['cpg_retention_readpos'].items() if 'no_data_available' not in dd.keys()]),
              dict([(sid, dd['2']) for sid, dd in self.mdata['cpg_retention_readpos'].items() if 'no_data_available' not in dd.keys()]),
              dict([(sid, dd['1']) for sid, dd in self.mdata['cph_retention_readpos'].items() if 'no_data_available' not in dd.keys()]),
              dict([(sid, dd['2']) for sid, dd in self.mdata['cph_retention_readpos'].items() if 'no_data_available' not in dd.keys()])]

        pconfig = {
            'id': 'biscuit_retention_cytosine',
            'title': 'BISCUIT: Retention vs. Base Position in Read',
            'xlab': 'Position in Read',
            'ylab': 'CpG Retention Rate (%)',
            'ymin': 0,
            'ymax': 100,
            'yMinRange': 0,
            'yFloor': 0,
            'tt_label': '<strong>Position {point.x}:</strong> {point.y:.2f}%',
            'data_labels': [{'name': 'CpG Read 1', 'ylab': 'CpG Retention Rate (%)'},
                            {'name': 'CpG Read 2', 'ylab': 'CpG Retention Rate (%)'},
                            {'name': 'CpH Read 1', 'ylab': 'CpH Retention Rate (%)'},
                            {'name': 'CpH Read 2', 'ylab': 'CpH Retention Rate (%)'}]
        }

        self.add_section(
            name = 'Retention vs. Base Position in Read',
            anchor = 'biscuit-retention-cytosine',
            description = 'Shows the distribution of cytosine retention rates ' \
            'across base positions in the read (a.k.a. the M-bias plot). Shown ' \
            'for cytosines in both a CpG and CpH context.',
            plot = linegraph.plot(pd, pconfig)
        )

    def parse_logs_cph_retention_readpos(self, f, fn):
        ''' Handled by parse_logs_cpg_retention_readpos() '''
        return self.parse_logs_cpg_retention_readpos(f, fn)

    def chart_cph_retention_readpos(self):
        ''' Handled by chart_cpg_retention_readpos() '''
        pass

    def parse_logs_read_avg_retention_rate(self, f, fn):
        '''
        Parses _totalReadConversionRate.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of read averaged fraction of retainied cytosines by context
        '''

        data = {}
        try:
            m = re.match(r'([\d.]+)\t([\d.]+)\t([\d.]+)\t([\d.]+)', f.splitlines()[2])
        except IndexError:
            return {}
        else:
            if m is not None:
                data['rca'] = 100.0 * float(m.group(1))
                data['rcc'] = 100.0 * float(m.group(2))
                data['rcg'] = 100.0 * float(m.group(3))
                data['rct'] = 100.0 * float(m.group(4))

        return data

    def chart_read_avg_retention_rate(self):
        '''
        Charts _totalReadConversionRate.txt
               _totalBaseConversionRate.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Retenion vs. Base Position in Read chart
        '''

        mdata_byread = {}
        for sid, dd in self.mdata['read_avg_retention_rate'].items():
            mdata_byread[sid] = dd

        mdata_bybase = {}
        for sid, dd in self.mdata['base_avg_retention_rate'].items():
            mdata_bybase[sid] = dd

        pdata = {}
        for sid, dd in mdata_byread.items():
            pdata[sid] = dict(list(dd.items()) + list(mdata_bybase[sid].items()))

        shared = {
            'format':'{:,.2f}',
            'min':0,
            'max':100,
            'suffix':'%',
            'scale': 'YlGnBu'
        }

        pheader = OrderedDict()
        pheader['rca'] = dict(shared, **{'title':'RA CpA', 'description':'Read Averaged CpA Retention'})
        pheader['rcc'] = dict(shared, **{'title':'RA CpC', 'description':'Read Averaged CpC Retention'})
        pheader['rcg'] = dict(shared, **{'title':'RA CpG', 'description':'Read Averaged CpG Retention'})
        pheader['rct'] = dict(shared, **{'title':'RA CpT', 'description':'Read Averaged CpT Retention'})
        pheader['bca'] = dict(shared, **{'title':'BA CpA', 'description':'Base Averaged CpA Retention'})
        pheader['bcc'] = dict(shared, **{'title':'BA CpC', 'description':'Base Averaged CpC Retention'})
        pheader['bcg'] = dict(shared, **{'title':'BA CpG', 'description':'Base Averaged CpG Retention'})
        pheader['bct'] = dict(shared, **{'title':'BA CpT', 'description':'Base Averaged CpT Retention'})

        pconfig = {
            'id': 'biscuit_retention',
            'table_title': 'BISCUIT: Cytosine Retention',
            'save_file': True,
            'sortRows': False
        }

        self.add_section(
            name = 'Cytosine Retention',
            anchor = 'biscuit-retention',
            description = 'Shows the cytosine retention rate for different ' \
            'contexts. `RA` stands for read-averaged rates and `BA` stands ' \
            'for base-averaged rates. See Help for more details.',
            helptext = 'Note, the cytosine retention rate is ' \
            '`1 - (cytosine conversion rate)` Additionally, assuming complete, ' \
            'but not over, bisulfite conversion, the cytosine retention rate ' \
            'is the average cytosine modification (including 5mC, 5hmC, ' \
            'etc) rate.',
            plot = table.plot(pdata, pheader, pconfig)
        )

    def parse_logs_base_avg_retention_rate(self, f, fn):
        '''
        Parses _totalBaseConversionRate.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of base averaged fraction of retainied cytosines by context
        '''

        data = {}
        try:
            m = re.match(r'([\d.]+)\t([\d.]+)\t([\d.]+)\t([\d.]+)', f.splitlines()[2])
        except IndexError:
            return {}
        else:
            if m is not None:
                data['bca'] = 100.0 * float(m.group(1))
                data['bcc'] = 100.0 * float(m.group(2))
                data['bcg'] = 100.0 * float(m.group(3))
                data['bct'] = 100.0 * float(m.group(4))

        return data

    def chart_base_avg_retention_rate(self):
        ''' Handled by chart_read_avg_retention_rate() '''
        pass
