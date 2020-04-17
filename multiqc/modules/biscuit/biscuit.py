#!/usr/bin/env python

''' MultiQC module to parse output from BISCUIT '''

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc.plots import linegraph, bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    '''
    Inherit from base module class
    Initializie data structures and prepare BISCUIT report
    Inputs:
        No inputs
    Returns:
        BISCUIT report for MultiQC
    '''
    def __init__(self):
        #Initialize the parent object
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
            'markdup_report': {},
            # Base coverage
            'covdist': {},
            'covdist_q40': {},
            'covdist_q40_botgc': {},
            'covdist_q40_topgc': {},
            # CpG coverage
            'covdist_cpg': {},
            'covdist_cpg_q40': {},
            'covdist_cpg_q40_botgc': {},
            'covdist_cpg_q40_topgc': {},
            # Uniformity
            'qc_cv': {},
            'qc_cpg_cv': {},
            # CpG distribution
            'qc_cpg_dist': {},
            # Retention
            'retention_dist': {},
            'retention_dist_byread': {},
            'retention_cpg_readpos': {},
            'retention_cph_readpos': {},
            'retention_rate_bybase': {},
            'retention_rate_byread': {}
        }

        file_suffixes = [
            # General statistics
            '_mapq_table',
            '_strand_table',
            '_isize_score_table',
            # Duplicate reporting
            '_dup_report',
            '_markdup_report',
            # Base coverage
            '_covdist_table',
            '_covdist_q40_table',
            '_covdist_q40_botgc_table',
            '_covdist_q40_topgc_table',
            # CpG coverage
            '_covdist_cpg_table',
            '_covdist_cpg_q40_table',
            '_covdist_cpg_q40_botgc_table',
            '_covdist_cpg_q40_topgc_table',
            # Uniformity
            '_all_cv_table',
            '_cpg_cv_table',
            # CpG distribution
            '_cpg_dist_table',
            # Retention
            '_CpGRetentionDist',
            '_freqOfTotalRetentionPerRead',
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

                # Clean uncommon file suffixes
                for file_suffix in file_suffixes:
                    s_name = s_name.replace(file_suffix, '')

                if s_name in self.mdata[k]:
                    log.debug('Duplicate sample name found in {}! Overwriting: {}'.format(f['fn'], s_name))

                self.mdata[k][s_name] = getattr(self, 'parse_logs_{}'.format(k))(f['f'], f['fn'])

        for k in self.mdata:
            self.mdata[k] = self.ignore_samples(self.mdata[k])

        if sum([len(self.mdata[k]) for k in self.mdata]) == 0:
            raise UserWarning

        # Basic stats table
        self.biscuit_stats_table()

        # Write out BISCUIT MultiQC report
        log.info('Found {} reports'.format(sum([len(self.mdata[k]) for k in self.mdata])))
        for k in self.mdata:
            if len(self.mdata[k]) > 0:
                self.write_data_file(self.mdata[k], 'multiqc_biscuit_{}'.format(k), data_format='json')
                log.debug('Found {} {} reports'.format(len(self.mdata[k]), k))
                getattr(self, 'chart_{}'.format(k))()

    def biscuit_stats_table(self):
        '''
        Create general statistics table for BISCUIT data
        Inputs:
            No true inputs, but
                needs mdata['align_mapq'] and mdata['dup_report'] filled
        Returns:
            No true returns, but creates stats table
        '''
        pd = {}

        for sid, dd in self.mdata['align_mapq'].items():
            if 'no_data_available' not in dd.keys():
                allreads = sum([int(_) for _ in dd.values()])
                pd[sid] = {'%aligned': 100.0 * float(allreads - int(dd['unmapped'])) / allreads}

        for sid, dd in self.mdata['dup_report'].items():
            if sid not in pd:
                pd[sid] = {}
            if 'all' in dd and dd['all'] != -1:
                pd[sid]['%dupRate_All'] = dd['all']
            if 'all-q40' in dd and dd['all'] != -1:
                pd[sid]['%dupRate_All-q40'] = dd['all-q40']

        pheader = OrderedDict()
        pheader['%aligned']         = {'title':'% Aligned', 'max':100, 'min':0, 'suffix':'%', 'scale':'Reds'}
        pheader['%dupRate_All']     = {'title':'% Overall Dup Rate', 'max':100, 'min':0, 'suffix':'%', 'scale':'Reds'}
        pheader['%dupRate_All-q40'] = {'title':'% Q40 Overall Dup Rate', 'max':100, 'min':0, 'suffix':'%', 'scale':'Purples'}
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
            return {'no_data_available': 1}

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
            No returns, generates Mapping Summary and Mapping Quality Distribution charts
        '''

        # Number of optimally and suboptimally mapped reads and unmapping reads
        pd = {}
        for sid, dd in self.mdata['align_mapq'].items():
            if 'no_data_available' not in dd.keys():
                pd[sid] = {'OAligned': 0, 'SAligned': 0, 'UAligned': 0}
                for mapq, cnt in dd.items():
                    if mapq == 'unmapped':
                        pd[sid]['UAligned'] += int(cnt)
                    elif int(mapq) >= 40:
                        pd[sid]['OAligned'] += int(cnt)
                    else:
                        pd[sid]['SAligned'] += int(cnt)

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
            plot = bargraph.plot(pd,
                                 OrderedDict([('OAligned', {'color': '#a6cee3', 'name': 'Optimally Aligned Reads'}),
                                              ('SAligned', {'color': '#1f78b4', 'name': 'Suboptimally Aligned Reads'}),
                                              ('UAligned', {'color': '#b2df8a', 'name': 'Unaligned Reads'})]),
                                 {'id': 'biscuit_mapping_overview',
                                  'title': 'BISCUIT: Mapping Overview',
                                  'ylab': 'Number of Reads',
                                  'cpswitch_counts_label': '# Reads'})
        )

        # Mapping quality distribution
        total = {}
        for sid, dd in self.mdata['align_mapq'].items():
            if 'no_data_available' not in dd.keys():
                total[sid] = sum([int(cnt) for _, cnt in dd.items() if _ != 'unmapped'])

        pd_mapping = {}
        for sid, dd in self.mdata['align_mapq'].items():
            if 'no_data_available' not in dd.keys():
                mapqcnts = []
                for mapq in range(61):
                    if str(mapq) in dd:
                        mapqcnts.append(100.0 * float(dd[str(mapq)]) / total[sid])
                    else:
                        mapqcnts.append(0)
                pd_mapping[sid] = dict(zip(range(61), mapqcnts))

        self.add_section(
            name = 'Mapping Quality Distribution',
            anchor = 'biscuit-mapq',
            description = 'Shows the percentage of the total number of mapped ' \
            'reads each mapping quality has (for primary alignments only).',
            plot = linegraph.plot(pd_mapping,
                                  {'id': 'biscuit_mapq',
                                   'title': 'BISCUIT: Distribution of Mapping Qualities',
                                   'ymin': 0,
                                   'xmin': 0,
                                   'yLabelFormat': '{value}%',
                                   'tt_label': '<strong>Q{point.x}:</strong> {point.y:.2f}% of reads',
                                   'ylab': '% of Primary Mapped Reads',
                                   'xlab': 'Mapping Quality'})
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
        file_data = f.splitlines()[1:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug('No data available in {}. Will not fill corresponding entries.'.format(fn))
            return {'no_data_available': 1}

        data = {'read1': {}, 'read2': {}}
        for l in file_data:
            m = re.search(r'strand\t([12])([+-]*)\t(\d+)', l)
            if m is not None:
                if m.group(1) == '1':
                    data['read1'][m.group(2)] = int(m.group(3))
                if m.group(1) == '2':
                    data['read2'][m.group(2)] = int(m.group(3))

        return data

    def chart_align_strand(self):
        '''
        Chart _strand_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Mapping Strand Distribution chart
        '''

        # Mapping strand distribution
        pd1 = {}
        pd2 = {}
        for sid, dd in self.mdata['align_strand'].items():
            if 'no_data_available' not in dd.keys():
                pd1[sid] = dd['read1']
                pd2[sid] = dd['read2']

        # TODO: When PBAT mode is implemented, add comment in help text about
        #       how to interpret PBAT mode results
        self.add_section(
            name = 'Mapping Strand Distribution',
            anchor = 'biscuit-strands',
            description = 'For primary alignments, shows the number of reads ' \
            'mapped to each strand. See Help for more details',
            helptext = 'Most bisulfite libraries map read 1 to the parent ' \
            'strand (`++` or `--`) and read 2 to the daughter/synthesized ' \
            'strand (`+-` or `-+`). PBAT and most single-cell/low input ' \
            'libraries often do not follow this assumption.',
            plot = bargraph.plot([pd1, pd2],
                                 [OrderedDict([('++', {'color': '#F53855', 'name': '++: Waston-Aligned, Waston-Bisulfite Conversion'}),
                                               ('+-', {'color': '#E37B40', 'name': '+-: Waston-Aligned, Crick-Bisulfite Conversion' }),
                                               ('-+', {'color': '#46B29D', 'name': '-+: Crick-Aligned, Waston-Bisulfite Conversion' }),
                                               ('--', {'color': '#324D5C', 'name': '--: Crick-Aligned, Crick-Bisulfite Conversion'  })]),
                                  OrderedDict([('++', {'color': '#F53855', 'name': '++: Waston-Aligned, Waston-Bisulfite Conversion'}),
                                               ('+-', {'color': '#E37B40', 'name': '+-: Waston-Aligned, Crick-Bisulfite Conversion' }),
                                               ('-+', {'color': '#46B29D', 'name': '-+: Crick-Aligned, Waston-Bisulfite Conversion' }),
                                               ('--', {'color': '#324D5C', 'name': '--: Crick-Aligned, Crick-Bisulfite Conversion'  })])],
                                 {'id': 'biscuit_strands',
                                  'title': 'BISCUIT: Mapping Strand Distribution',
                                  'ylab': 'Number of Reads',
                                  'cpswitch_counts_label': '# Reads',
                                  'data_labels': [{'name': 'Read 1'},
                                                  {'name': 'Read 2'}]})
        )

    def parse_logs_align_isize(self, f, fn):
        '''
        Parse _isize_score_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of insert size and SW-score data
        '''
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug('No data available in {}. Will not fill corresponding entries.'.format(fn))
            return {'no_data_available': 1}

        data = {'I': {}, 'S': {}} # I = insert size, S = SW-score
        for l in file_data:
            fields = l.split()
            if fields[0] == 'I':
                data[fields[0]][int(fields[1])] = float(fields[2]) * 100.0
            elif fields[0] == 'S':
                data[fields[0]][fields[1]] = float(fields[2])

        return data

    def chart_align_isize(self):
        '''
        Chart _isize_score_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Insert Size Distribution chart
        '''

        pd_isize = {}
        for sid, dd in self.mdata['align_isize'].items():
            if 'no_data_available' not in dd.keys():
                pd_isize[sid] = dd['I']

        self.add_section(
            name = 'Insert Size Distribution',
            anchor = 'biscuit-isize',
            description = 'Shows the distribution of insert sizes. See Help ' \
            'for more details.',
            helptext = 'Insert size is defined as: `(right-most coordinate ' \
            'of reverse-mate read) - (left-most coordinate of forward-mate ' \
            'read)`. Insert sizes are calculated for reads with a "mapped in ' \
            'proper pair" `samtools` flag, an alignment score >= 40, and ' \
            'MAPQ>=40.',
            plot = linegraph.plot(pd_isize,
                                  {'id': 'biscuit_isize',
                                   'title': 'BISCUIT: Insert Size Distribution',
                                   'ymin': 0,
                                   'xmin': 0,
                                   'yLabelFormat': '{value}%',
                                   'smooth_points': 500, # limit number of points / smooth data
                                   'tt_label': '<strong>IS{point.x}:</strong> {point.y:.2f}% of reads',
                                   'ylab': '% Mapped Reads',
                                   'xlab': 'Insert Size'})
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
            (r'#bases covered by all reads: (\d+)',
             r'#bases covered by duplicate reads: (\d+)', 'all'),
            (r'#high-GC bases covered by all reads: (\d+)',
             r'#high-GC bases covered by duplicate reads: (\d+)', 'topGC'),
            (r'#low-GC bases covered by all reads: (\d+)',
             r'#low-GC bases covered by duplicate reads: (\d+)', 'lowGC'),
            (r'#bases covered by all q40-reads: (\d+)',
             r'#bases covered by duplicate q40-reads: (\d+)', 'all-q40'),
            (r'#high-GC bases covered by all q40-reads: (\d+)',
             r'#high-GC bases covered by duplicate q40-reads: (\d+)', 'topGC-q40'),
            (r'#low-GC bases covered by all q40-reads: (\d+)',
             r'#low-GC bases covered by duplicate q40-reads: (\d+)', 'botGC-q40')]

        data = {}
        for pat_all, pat_dup, k in patterns:
            m1 = re.search(pat_all, f, re.MULTILINE)
            m2 = re.search(pat_dup, f, re.MULTILINE)
            if m1 is not None and m2 is not None and float(m1.group(1)) > 0:
                data[k] = (100.0 * float(m2.group(1)) / float(m1.group(1)))
            else:
                data[k] = -1

        return data

    def chart_dup_report(self):
        '''
        Charts _dup_report.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates CpG Coverage by Genomic Feature
            and CpG Island Coverage charts
        '''

        pd = dict([(sid, dd) for sid, dd in self.mdata['dup_report'].items() if dd['all'] != -1])

        self.add_section(
            name = 'Duplicate Rates',
            anchor = 'biscuit-dup-report',
            description = 'Shows the percentage of bases that are duplicates ' \
            'out of the total number of bases. See Help for more details.',
            helptext = '"High GC" and "low GC" content regions are the top ' \
            'and bottom 10% of 100bp windows for GC content, respectively.',
            plot = table.plot(pd,
                              OrderedDict([('all', {'title': 'Overall', 'suffix': '%', 'max': 100, 'min': 0, 'scale': 'Reds'}),
                                           ('topGC', {'title': 'Overall High GC', 'suffix': '%', 'max': 100, 'min': 0, 'scale': 'PuRd'}),
                                           ('lowGC', {'title': 'Overall Low GC', 'suffix': '%', 'max': 100, 'min': 0, 'scale': 'PuRd'}),
                                           ('all-q40', {'title': 'Q40', 'suffix': '%', 'max': 100, 'min': 0, 'scale': 'Purples'}),
                                           ('topGC-q40', {'title': 'Q40 High GC', 'suffix': '%', 'max': 100, 'min': 0, 'scale': 'PuBuGn'}),
                                           ('botGC-q40', {'title': 'Q40 Low GC', 'suffix': '%', 'max': 100, 'min': 0, 'scale': 'PuBuGn'})]),
                              {'id': 'biscuit_dup_report',
                               'table_title': 'BISCUIT: Duplicate Rates',
                               'save_file': True,
                               'sortRows': False})
        )

    # TODO: Remove when biscuit markdup is officially removed from code base
    def parse_logs_markdup_report(self, f, fn):
        '''
        Parses _markdup_report.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of duplicate rates from biscuit markdup
        '''
        data = {}
        m = re.search(r'marked (\d+) duplicates from (\d+) paired-end reads', f, re.MULTILINE)
        if m is None:
            data = {'PE': None, 'PEdup': None, 'dupRatePE': None}
        else:
            data['PEdup'] = int(m.group(1))
            data['PE'] = int(m.group(2))
            if data['PE'] > 0:
                data['dupRatePE'] = 100.0 * float(data['PEdup']) / data['PE']
            else:
                data['dupRatePE'] = None

        m = re.search(r'marked (\d+) duplicates from (\d+) single-end reads', f, re.MULTILINE)
        if m is None:
            data = {'SE': None, 'SEdup': None, 'dupRateSE': None}
        else:
            data['SEdup'] = int(m.group(1))
            data['SE'] = int(m.group(2))
            if data['SE'] > 0:
                data['dupRateSE'] = 100.0 * float(data['SEdup']) / data['SE']
            else:
                data['dupRateSE'] = None
                
        m = re.search(r'identified (\d+) dangling paired-end reads', f, re.MULTILINE)
        if m is None:
            data['PEdangling'] = None
        else:
            data['PEdangling'] = int(m.group(1))
        
        return data

    def chart_markdup_report(self):
        ''' biscuit markdup is being deprecated, so pass until officially removed '''
        pass

    ########################################
    #### Base Coverage and CpG Coverage ####
    ########################################
    def parse_logs_covdist(self, f, fn):
        '''
        Parses _covdist_table.txt
               _covdist_q40_table.txt
               _covdist_q40_botgc_table.txt
               _covdist_q40_topgc_table.txt
               _covdist_cpg_table.txt
               _covdist_cpg_q40_table.txt
               _covdist_cpg_q40_botgc_table.txt
               _covdist_cpg_q40_topgc_table.txt
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
            ccov_cnts.append(_ccov_cnt/1000000.)
            _ccov_cnt -= dd[cov]

        return dict(zip(covs, ccov_cnts))

    def parse_logs_covdist_q40(self, f, fn):
        ''' Handled by parse_logs_covdist() '''
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_q40_botgc(self, f, fn):
        ''' Handled by parse_logs_covdist() '''
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_q40_topgc(self, f, fn):
        ''' Handled by parse_logs_covdist() '''
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_cpg(self, f, fn):
        ''' Handled by parse_logs_covdist() '''
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_cpg_q40(self, f, fn):
        ''' Handled by parse_logs_covdist() '''
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_cpg_q40_botgc(self, f, fn):
        ''' Handled by parse_logs_covdist() '''
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_cpg_q40_topgc(self, f, fn):
        ''' Handled by parse_logs_covdist() '''
        return self.parse_logs_covdist(f, fn)

    def chart_covdist(self):
        '''
        Charts _covdist_table.txt
               _covdist_q40_table.txt
               _covdist_q40_botgc_table.txt
               _covdist_q40_topgc_table.txt
               _covdist_cpg_table.txt
               _covdist_cpg_q40_table.txt
               _covdist_cpg_q40_botgc_table.txt
               _covdist_cpg_q40_topgc_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Overall Coverage chart
        '''

        mdata = [
            self.mdata['covdist'],
            self.mdata['covdist_q40'],
            self.mdata['covdist_q40_botgc'],
            self.mdata['covdist_q40_topgc'],
            self.mdata['covdist_cpg'],
            self.mdata['covdist_cpg_q40'],
            self.mdata['covdist_cpg_q40_botgc'],
            self.mdata['covdist_cpg_q40_topgc']
        ]

        self.add_section(
            name = 'Cumulative Coverage',
            anchor = 'biscuit-cumulative-coverage',
            description = 'Shows the number of reads covering all bases or ' \
            'CpGs only. See Help for more details.',
            helptext = '"Q40" shows only those reads with mapping quality ' \
            'greater than or equal to 40. "High GC" and "low GC" content ' \
            'regions are the top and bottom 10% of 100bp windows for GC ' \
            'content, respectively.',
            plot = linegraph.plot(mdata,
                                  {'id': 'biscuit_cumulative',
                                   'title': 'BISCUIT: Cumulative Coverage',
                                   'ymin': 0,
                                   'xLabelFormat': '{value}X',
                                   'tt_label': '<strong>{point.x}X:</strong> {point.y:.2f}M',
                                   'xlab': 'Sequencing Depth',
                                   'data_labels': [{'name': 'Bases (all)', 'ylab':'Millions of Bases'},
                                                   {'name': 'Bases Q40', 'ylab':'Millions of Bases (Q40)'},
                                                   {'name': 'Bases Q40 low GC', 'ylab':'Millions of Low-GC Bases (Q40)'},
                                                   {'name': 'Bases Q40 high GC', 'ylab':'Millions of High-GC Bases (Q40)'},
                                                   {'name': 'CpG (all)', 'ylab':'Millions of CpGs'},
                                                   {'name': 'CpG Q40', 'ylab':'Millions of CpGs (Q40)'},
                                                   {'name': 'CpG Q40 low GC', 'ylab':'Millions of Low-GC CpGs (Q40)'},
                                                   {'name': 'CpG Q40 high GC', 'ylab':'Millions of High-GC CpGs (Q40)'}]})
        )

        basecov = OrderedDict()
        for sid, dd in mdata[0].items():
            if sid not in basecov:
                basecov[sid] = {}

            if dd[0] != -1:
                basecov[sid]['all'] = float(dd[1])/float(dd[0])*100.0

        for sid, dd in mdata[1].items():
            if sid not in basecov:
                basecov[sid] = {}

            if dd[0] != -1:
                basecov[sid]['q40'] = float(dd[1])/float(dd[0])*100.0

        for sid, dd in mdata[4].items():
            if sid not in basecov:
                basecov[sid] = {}

            if dd[0] != -1:
                basecov[sid]['cpg_all'] = float(dd[1])/float(dd[0])*100.0

        for sid, dd in mdata[5].items():
            if sid not in basecov:
                basecov[sid] = {}

            if dd[0] != -1:
                basecov[sid]['cpg_q40'] = float(dd[1])/float(dd[0])*100.0

        # base coverage >=1x table
        if len(basecov) == 0:
            return

        self.add_section(
            name = 'Coverage by at Least One Read',
            anchor = 'biscuit-one-read-coverage',
            description = 'The fraction of the genome and genomic CpGs ' \
            'covered by at least one read. See Help for more details',
            helptext = '"Q40" shows the percentage for reads with MAPQ>=40',
            plot = table.plot(basecov,
                              OrderedDict([('all', {'title': 'Genome (All)', 'max': 100, 'min': 0, 'suffix': '%', 'scale': 'Reds'}),
                                           ('q40', {'title': 'Genome (Q40)', 'max': 100, 'min': 0, 'suffix': '%', 'scale': 'Purples'}),
                                           ('cpg_all', {'title': 'CpG (All)', 'max': 100, 'min': 0, 'suffix': '%', 'scale': 'Blues'}),
                                           ('cpg_q40', {'title': 'CpG (Q40)', 'max': 100, 'min': 0, 'suffix': '%', 'scale': 'Oranges'})]),
                              {'id': 'biscuit_one_read_coverage',
                               'table_title': 'BISCUIT: Coverage by at Least One Read',
                               'save_file': True,
                               'sortRows': False})
        )

    def chart_covdist_q40(self):
        ''' Handled in chart_covdist() '''
        pass

    def chart_covdist_q40_botgc(self):
        ''' Handled in chart_covdist() '''
        pass

    def chart_covdist_q40_topgc(self):
        ''' Handled in chart_covdist() '''
        pass

    def chart_covdist_cpg(self):
        ''' Handled in chart_covdist() '''
        pass

    def chart_covdist_cpg_q40(self):
        ''' Handled in chart_covdist() '''
        pass

    def chart_covdist_cpg_q40_botgc(self):
        ''' Handled in chart_covdist() '''
        pass

    def chart_covdist_cpg_q40_topgc(self):
        ''' Handled in chart_covdist() '''
        pass

    ########################################
    ####      Depths and Uniformity     ####
    ########################################
    def parse_logs_qc_cv(self, f, fn):
        '''
        Parses _all_cv_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of depth uniformity measures
        '''

        data = {}
        targets = ['all', 'all_topgc', 'all_botgc']
        for k in targets:
            m = re.search('_{}\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)'.format(k),
                          f, re.MULTILINE)
            if m is None:
                data[k] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data[k] = {'mu': float(m.group(1)),
                           'sigma': float(m.group(2)),
                           'cv': float(m.group(3))}

        return data

    def chart_qc_cv(self):
        '''
        Charts _all_cv_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Sequencing Depth - Whole Genome chart
        '''

        # Whole Genome - Sequencing depth and uniformity
        pd_wg = OrderedDict()
        for sid, dd in self.mdata['qc_cv'].items():
            pd_wg[sid] = OrderedDict()
            for ctg in ['all', 'all_topgc', 'all_botgc']:
                if ctg in dd:
                    if dd[ctg]['mu'] != -1:
                        pd_wg[sid]['mu_'+ctg] = dd[ctg]['mu']
                    if dd[ctg]['cv'] != -1:
                        pd_wg[sid]['cv_'+ctg] = dd[ctg]['cv']

        if len(pd_wg) == 0:
            return

        pheader_wg = OrderedDict()
        pheader_wg['mu_all'] = {'title': 'Genome Mean', 'description': 'Whole Genome Mean', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'Reds'}
        pheader_wg['mu_all_topgc'] = {'title': 'High GC Mean', 'description': 'Top Decile for GC Content Mean', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'PuRd'}
        pheader_wg['mu_all_botgc'] = {'title': 'Low GC Mean', 'description': 'Bottom Decile for GC Content Mean', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'PuRd'}
        pheader_wg['cv_all'] = {'title': 'Genome CV', 'description': 'Whole Genome Coeff. of Var.', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'Reds'}
        pheader_wg['cv_all_topgc'] = {'title': 'High GC CV', 'description': 'Top Decile for GC Content Coeff. of Var.', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'PuRd'}
        pheader_wg['cv_all_botgc'] = {'title': 'Low GC CV', 'description': 'Bottom Decile for GC Content Coeff. of Var.', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'PuRd'}

        self.add_section(
            name = 'Sequencing Depth - Whole Genome',
            anchor = 'biscuit-seq-depth-all',
            description = 'Shows the sequence depth mean and uniformity ' \
            'measured by the Coefficient of Variation (CV, defined as ' \
            'std. dev./mean). See Help for more details.',
            helptext = 'Mean and CV are calculated for reads with MAPQ>=40 ' \
            'across the entire genome. "High GC" and "Low GC" are defined as ' \
            'the top and bottom deciles for GC content as measured in 100bp ' \
            'non-overlapping windows, respectively.',
            plot = table.plot(pd_wg,
                              pheader_wg,
                              {'id': 'biscuit_seq_depth_all',
                               'table_title': 'BISCUIT: Sequencing Depth - Whole Genome',
                               'save_file': True,
                               'sortRows': False})
        )

    def parse_logs_qc_cpg_cv(self, f, fn):
        '''
        Parses _cpg_cv_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of depth uniformity measures
        '''

        data = {}
        targets = ['cpg', 'cpg_topgc', 'cpg_botgc']
        for k in targets:
            m = re.search('_{}\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)'.format(k),
                          f, re.MULTILINE)
            if m is None:
                data[k] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data[k] = {'mu': float(m.group(1)),
                           'sigma': float(m.group(2)),
                           'cv': float(m.group(3))}

        return data

    def chart_qc_cpg_cv(self):
        '''
        Charts _cpg_cv_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Sequencing Depth - CpGs Only chart
        '''

        # CpGs Only - Sequencing depth and uniformity
        pd_cg = OrderedDict()
        for sid, dd in self.mdata['qc_cpg_cv'].items():
            pd_cg[sid] = OrderedDict()
            for ctg in ['cpg', 'cpg_topgc', 'cpg_botgc']:
                if ctg in dd:
                    if dd[ctg]['mu'] != -1:
                        pd_cg[sid]['mu_'+ctg] = dd[ctg]['mu']
                    if dd[ctg]['cv'] != -1:
                        pd_cg[sid]['cv_'+ctg] = dd[ctg]['cv']

        if len(pd_cg) == 0:
            return

        pheader_cg = OrderedDict()
        pheader_cg['mu_cpg'] = {'title': 'All CpG Mean', 'description': 'All CpGs Mean', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'Oranges'}
        pheader_cg['mu_cpg_topgc'] = {'title': 'High GC CpG Mean', 'description': 'CpGs in Top Decile for GC Content Mean', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'YlOrRd'}
        pheader_cg['mu_cpg_botgc'] = {'title': 'Low GC CpG Mean', 'description': 'CpGs in Bottom Decile for GC Content Mean', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'YlOrRd'}
        pheader_cg['cv_cpg'] = {'title': 'All CpG CV', 'description': 'All CpGs Coeff. of Var.', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'Oranges'}
        pheader_cg['cv_cpg_topgc'] = {'title': 'High GC CpG CV', 'description': 'CpGs in Top Decile for GC Content Coeff. of Var.', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'YlOrRd'}
        pheader_cg['cv_cpg_botgc'] = {'title': 'Low GC CpG CV', 'description': 'CpGs in Bottom Decile for GC Content Coeff. of Var.', 'max': 50, 'min': 0, 'format': '{:,.3f}', 'scale': 'YlOrRd'}

        self.add_section(
            name = 'Sequencing Depth - CpGs Only',
            anchor = 'biscuit-seq-depth-cpg',
            description = 'Shows the sequence depth mean and uniformity ' \
            'measured by the Coefficient of Variation (CV, defined as ' \
            'std. dev./mean) for CpG locations. See Help for more details.',
            helptext = 'Mean and CV are calculated for CpGs in reads with ' \
            'MAPQ>=40 across the entire genome. "High GC" and "Low GC" are ' \
            'defined as the top and bottom deciles for GC content as measured ' \
            'in 100bp non-overlapping windows, respectively.',
            plot = table.plot(pd_cg,
                              pheader_cg,
                              {'id': 'biscuit_seq_depth_cpg',
                               'table_title': 'BISCUIT: Sequencing Depth - CpGs Only',
                               'save_file': True,
                               'sortRows': False})
        )

    ########################################
    ####        CpG Distribution        ####
    ########################################
    def parse_logs_qc_cpg_dist(self, f, fn):
        '''
        Parses _cpg_dist_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of CpG coverage distributions
        '''

        data = {}
        for ctg in ['TotalCpGs', 'ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
            m = re.search('{}\t(\d+)\t(\d+)\t(\d+)'.format(ctg), f, re.MULTILINE)
            if m is None:
                data[ctg] = {'a': -1, 'uc': -1, 'ac': -1}
            else:
                data[ctg] = {'a': int(m.group(1)),
                             'uc': int(m.group(2)),
                             'ac': int(m.group(3))}

        num_cgi = -999
        m = re.search(r'#CpG Islands\t(\d+)', f, re.MULTILINE)
        if m is None:
            num_cgi = -1
        else:
            num_cgi = int(m.group(1))

        patterns = [
            (r'one CpG\t(\d+)', 'one'),
            (r'three CpGs\t(\d+)', 'three'),
            (r'five CpGs\t(\d+)', 'five'),
            (r'ten CpGs\t(\d+)', 'ten')]

        data['cgi_coverage'] = {}
        for pat, k in patterns:
            m = re.search(pat, f, re.MULTILINE)
            if m is None or num_cgi == -1:
                data['cgi_coverage'][k] = -1
            else:
                data['cgi_coverage'][k] = 100.0 * float(m.group(1)) / num_cgi

        return data

    def chart_qc_cpg_dist(self):
        '''
        Charts _cpg_dist_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates CpG Coverage by Genomic Feature
            and CpG Island Coverage charts
        '''

        if len(self.mdata['qc_cpg_dist']) == 0:
            return

        # Assorted regions
        pd = OrderedDict()

        # TODO: Try to figure out what is meant to be shown here
        # If figured out, put back in and update add_section description note
        #pd['Genome'] = OrderedDict()
        #dd = list(self.mdata['qc_cpg_dist'].values())[0]
        #for ctg in ['ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
        #    ctg1 = ctg.replace('CpGs','')
        #    pd['Genome'][ctg1] = 100.0 * float(dd[ctg]['uc']) / dd['TotalCpGs']['uc']

        for sid, dd in self.mdata['qc_cpg_dist'].items():
            pd[sid] = OrderedDict()
            for ctg in ['ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
                if dd[ctg]['uc'] != -1:
                    ctg1 = ctg.replace('CpGs','')
                    pd[sid][ctg1] = 100.0 * float(dd[ctg]['uc']) / dd['TotalCpGs']['uc']

        self.add_section(
            name = 'CpG Coverage by Genomic Feature',
            anchor = 'biscuit-coverage-cpg-dist',
            description = 'Shows the fraction of uniquely covered CpGs for ' \
            'different categories relative to the total number of uniquely ' \
            'covered CpGs in the dataset. CpGs are from reads with MAPQ>=40.',
            plot = table.plot(pd,
                              OrderedDict([('Exonic', {'title': 'Exonic CpGs', 'max': 100, 'min': 0, 'suffix': '%', 'scale': 'Oranges'}),
                                           ('Repeat', {'title': 'Repeat-Masked CpGs', 'max': 100, 'min': 0, 'suffix': '%', 'scale': 'Oranges'}),
                                           ('Genic', {'title': 'Genic CpGs', 'max': 100, 'min': 0, 'suffix': '%', 'scale': 'Oranges'}),
                                           ('CGI', {'title': 'CpG Island CpGs', 'max': 100, 'min': 0, 'suffix': '%', 'scale': 'Oranges'})]),
                              {'id': 'biscuit_coverage_cpg_dist',
                               'table_title': 'BISCUIT: CpG Coverage by Genomic Feature',
                               'save_file': True,
                               'sortRows': False})
        )

        # CpG Islands
        pd = dict([(sid, dd['cgi_coverage']) for sid, dd in self.mdata['qc_cpg_dist'].items() if dd['cgi_coverage']['one'] != -1])

        self.add_section(
            name = 'CpG Island Coverage',
            anchor = 'biscuit-coverage-cgi',
            description = 'Shows the percentage of CpG islands (out of all ' \
            'CpG islands in the genome) that have at least 1, 3, 5, or 10 ' \
            'different CpGs covered. Coverage is based on reads with MAPQ>=40.',
            plot = table.plot(pd,
                              OrderedDict([('one', {'title': '>=1', 'max': 100, 'min': 0, 'suffix': '%', 'description': 'CpG islands with at least one CpG covered', 'scale': 'Oranges'}),
                                           ('three', {'title': '>=3', 'max': 100, 'min': 0, 'suffix': '%', 'description': 'CpG islands with at least three CpGs covered', 'scale': 'Oranges'}),
                                           ('five', {'title': '>=5', 'max': 100, 'min': 0, 'suffix': '%', 'description': 'CpG islands with at least five CpGs covered', 'scale': 'Oranges'}),
                                           ('ten', {'title': '>=10', 'max': 100, 'min': 0, 'suffix': '%', 'description': 'CpG islands with at least ten CpGs covered', 'scale': 'Oranges'})]),
                              {'id': 'biscuit_coverage_cgi',
                               'table_title': 'BISCUIT: CpG Island Coverage',
                               'save_file': True,
                               'sortRows': False})
        )

    ########################################
    ####          CpG Retention         ####
    ########################################
    def parse_logs_retention_dist(self, f, fn):
        '''
        Parses _CpGRetentionDist.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of distribution of CpG retention
        '''
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug('No data available in {}. Will not fill corresponding entries.'.format(fn))
            return {'no_data_available': 1}

        sumcounts = 0
        for l in file_data:
            fields = l.split('\t')
            if fields[0] != -1:
                sumcounts += int(fields[1])

        data = {}
        for l in file_data:
            fields = l.split('\t')
            if fields[0] != -1:
                data[int(fields[0])] = int(fields[1]) / float(sumcounts)

        if len(data) == 0:
            log.debug('Only low coverage entries in {}. Will not fill corresponding entries.'.format(fn))
            return {'no_data_available': 1}

        return data

    def chart_retention_dist(self):
        '''
        Charts _CpGRetentionDist.txt
               _freqOfTotalRetentionPerRead.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Number of Retention Distribution chart
        '''

        mdata_meth = {}
        for sid, dd in self.mdata['retention_dist'].items():
            if 'no_data_available' not in dd.keys():
                mdata_meth[sid] = dd

        mdata = self.mdata['retention_dist_byread']

        pd = [mdata_meth,
              dict([(sid, dd['CA']) for sid, dd in mdata.items() if 'no_data_available' not in dd.keys()]),
              dict([(sid, dd['CC']) for sid, dd in mdata.items() if 'no_data_available' not in dd.keys()]),
              dict([(sid, dd['CG']) for sid, dd in mdata.items() if 'no_data_available' not in dd.keys()]),
              dict([(sid, dd['CT']) for sid, dd in mdata.items() if 'no_data_available' not in dd.keys()])]

        # TODO: Improve description on plot
        self.add_section(
            name = 'Number of Retention Distribution',
            anchor = 'biscuit-retention-read',
            description = 'The "CpG Retention" tab shows the cytosine retention ' \
            'percentage and the fraction of CpGs having this percentage. ' \
            'The "Within-read CpN" tabs show how many reads have a given ' \
            'number of retained cytosines per read. See Help for more details.',
            helptext = 'The "Within-read CpN" tabs are truncated at 10 retained ' \
            'cytosines for ease of viewing.',
            plot = linegraph.plot(pd,
                                  {'id': 'biscuit_retention_read',
                                   'title': 'BISCUIT: Retention Distribution',
                                   'xlab': 'Number of Retained Cytosines Within Read',
                                   'data_labels': [{'name': 'CpG Retention', 'ylab': 'Fraction of CpGs', 'xlab': 'Cytosine Retention Percentage'},
                                                   {'name': 'Within-read CpA', 'ylab': 'Millions of Reads'},
                                                   {'name': 'Within-read CpC', 'ylab': 'Millions of Reads'},
                                                   {'name': 'Within-read CpG', 'ylab': 'Millions of Reads'},
                                                   {'name': 'Within-read CpT', 'ylab': 'Millions of Reads'}]})
        )

    def parse_logs_retention_dist_byread(self, f, fn):
        '''
        Parses _freqOfTotalRetentionPerRead.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of distribution of CpG retention by read in a given context
        '''
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug('No data available in {}. Will not fill corresponding entries.'.format(fn))
            return {'no_data_available': 1}

        data = {'CA': {}, 'CC': {}, 'CG': {}, 'CT': {}}
        for l in file_data:
            fields = l.strip().split('\t')
            if fields[0] not in data:
                return {}
            if int(fields[1]) <= 10: # remove number retained greater than 10
                data[fields[0]][int(fields[1])] = float(fields[2]) / 1000000.

        return data

    def chart_retention_dist_byread(self):
        ''' Handled in chart_retention_dist() '''
        pass

    def parse_logs_retention_cpg_readpos(self, f, fn):
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

    def chart_retention_cpg_readpos(self):
        '''
        Charts _CpGRetentionByReadPos.txt
               _CpHRetentionByReadPos.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Retenion vs. Base Position in Read chart
        '''
        
        mdata = [dict([(sid, dd['1']) for sid, dd in self.mdata['retention_cpg_readpos'].items() if 'no_data_available' not in dd.keys()]),
                 dict([(sid, dd['2']) for sid, dd in self.mdata['retention_cpg_readpos'].items() if 'no_data_available' not in dd.keys()]),
                 dict([(sid, dd['1']) for sid, dd in self.mdata['retention_cph_readpos'].items() if 'no_data_available' not in dd.keys()]),
                 dict([(sid, dd['2']) for sid, dd in self.mdata['retention_cph_readpos'].items() if 'no_data_available' not in dd.keys()])]

        self.add_section(
            name = 'Retention vs. Base Position in Read',
            anchor = 'biscuit-retention-cytosine',
            description = 'Shows the distribution of cytosine retention rates ' \
            'across base positions in the read (a.k.a. the M-bias plot). Shown ' \
            'for cytosines in both a CpG and CpH context.',
            plot = linegraph.plot(mdata,
                                  {'id': 'biscuit_retention_cytosine',
                                   'title': 'BISCUIT: Retention vs. Base Position in Read',
                                   'xlab': 'Position in Read',
                                   'ymin': 0,
                                   'ymax': 100,
                                   'yMinRange': 0,
                                   'yFloor': 0,
                                   'tt_label': '<strong>Position {point.x}:</strong> {point.y:.2f}%',
                                   'data_labels': [{'name': 'CpG Read 1', 'ylab': 'CpG Retention Rate (%)'},
                                                   {'name': 'CpG Read 2', 'ylab': 'CpG Retention Rate (%)'},
                                                   {'name': 'CpH Read 1', 'ylab': 'CpH Retention Rate (%)'},
                                                   {'name': 'CpH Read 2', 'ylab': 'CpH Retention Rate (%)'}]})
        )

    def parse_logs_retention_cph_readpos(self, f, fn):
        ''' Handled by parse_logs_retention_cpg_readpos() '''
        return self.parse_logs_retention_cpg_readpos(f, fn)

    def chart_retention_cph_readpos(self):
        ''' Handled by chart_retention_cpg_readpos() '''
        pass

    def parse_logs_retention_rate_byread(self, f, fn):
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

    def chart_retention_rate_byread(self):
        '''
        Charts _totalReadConversionRate.txt
               _totalBaseConversionRate.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Retenion vs. Base Position in Read chart
        '''

        mdata_byread = {}
        for sid, dd in self.mdata['retention_rate_byread'].items():
            mdata_byread[sid] = dd

        mdata_bybase = {}
        for sid, dd in self.mdata['retention_rate_bybase'].items():
            mdata_bybase[sid] = dd

        pdata = {}
        for sid, dd in mdata_byread.items():
            pdata[sid] = dict(list(dd.items()) + list(mdata_bybase[sid].items()))

        pheader = OrderedDict()
        pheader['rca'] = {'title':'RA CpA', 'format':'{:,.2f}', 'min':0, 'max':100, 'description':'Read Averaged CpA Retention', 'suffix':'%', 'scale': 'YlGnBu'}
        pheader['rcc'] = {'title':'RA CpC', 'format':'{:,.2f}', 'min':0, 'max':100, 'description':'Read Averaged CpC Retention', 'suffix':'%', 'scale': 'YlGnBu'}
        pheader['rcg'] = {'title':'RA CpG', 'format':'{:,.2f}', 'min':0, 'max':100, 'description':'Read Averaged CpG Retention', 'suffix':'%', 'scale': 'YlGnBu'}
        pheader['rct'] = {'title':'RA CpT', 'format':'{:,.2f}', 'min':0, 'max':100, 'description':'Read Averaged CpT Retention', 'suffix':'%', 'scale': 'YlGnBu'}
        pheader['bca'] = {'title':'BA CpA', 'format':'{:,.2f}', 'min':0, 'max':100, 'description':'Base Averaged CpA Retention', 'suffix':'%', 'scale': 'YlGnBu'}
        pheader['bcc'] = {'title':'BA CpC', 'format':'{:,.2f}', 'min':0, 'max':100, 'description':'Base Averaged CpC Retention', 'suffix':'%', 'scale': 'YlGnBu'}
        pheader['bcg'] = {'title':'BA CpG', 'format':'{:,.2f}', 'min':0, 'max':100, 'description':'Base Averaged CpG Retention', 'suffix':'%', 'scale': 'YlGnBu'}
        pheader['bct'] = {'title':'BA CpT', 'format':'{:,.2f}', 'min':0, 'max':100, 'description':'Base Averaged CpT Retention', 'suffix':'%', 'scale': 'YlGnBu'}

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
            plot = table.plot(pdata,
                              pheader,
                              {'id': 'biscuit_retention',
                               'table_title': 'BISCUIT: Cytosine Retention',
                               'save_file': True,
                               'sortRows': False})
        )

    def parse_logs_retention_rate_bybase(self, f, fn):
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

    def chart_retention_rate_bybase(self):
        ''' Handled by chart_retention_rate_byread() '''
        pass
