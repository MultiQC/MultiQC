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
            info='is a tool to map bisulfite converted sequence reads and' \
            ' determine cytosine methylation states.')

        # Set up data structures
        self.mdata = {
            # General statistics
            'align_mapq': {},
            'align_strand': {},
            'align_isize': {},
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
            #'qc_cv': {},
            #'qc_cpg_cv': {},
            # CpG distribution
            #'qc_cpg_dist': {},
            # Duplicate reporting
            #'dup_report': {},
            #'markdup_report': {},
            # Retention
            #'retention_dist': {},
            #'retention_dist_byread': {},
            #'retention_cpg_readpos': {},
            #'retention_cph_readpos': {},
            #'retention_rate_bybase': {},
            #'retention_rate_byread': {}
        }

        file_suffixes = [
            # General statistics
            '_mapq_table',
            '_strand_table',
            '_isize_score_table',
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
            # Duplicate reporting
            '_dup_report',
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
                # Clean s_name before further processing
                s_name = self.clean_s_name(f['s_name'], f['root'])

                # Clean uncommon file suffixes
                for file_suffix in file_suffixes:
                    s_name = s_name.replace(file_suffix, '')

                if s_name in self.mdata[k]:
                    log.debug('Duplicate sample name found! Overwriting: {}'.format(s_name))

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
                self.write_data_file(self.mdata[k], 'multiqc_biscuit_{}'.format(k))
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
            if 'data_is_missing' not in dd.keys():
                allreads = sum([int(_) for _ in dd.values()])
                pd[sid] = {'%aligned': 100 * float(allreads - int(dd['unmapped'])) / allreads}

        #for sid, dd in self.mdata['dup_report'].items():
        #    if sid not in pd:
        #        pd[sid] = {}
        #    if 'all' in dd and dd['all'] is not None:
        #        pd[sid]['%dupRate_All'] = dd['all']
        #    if 'all-q40' in dd and dd['all-q40'] is not None:
        #        pd[sid]['%dupRate_All-q40'] = dd['all-q40']

        pheader = OrderedDict()
        pheader['%aligned']          = {'title':'% Aligned', 'max':100, 'min':0, 'suffix':'%', 'scale':'Greens'}
        #pheader['%dupRate_All']     = {'title':'% Overall Dup Rate', 'max':100, 'min':0, 'suffix':'%', 'scale':'Reds'}
        #pheader['%dupRate_All-q40'] = {'title':'% Q40 Overall Dup Rate', 'max':100, 'min':0, 'suffix':'%', 'scale':'Reds'}
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
            log.warning('Missing data in {}. Will not fill corresponding entries.'.format(fn))
            return {'data_is_missing': 1}

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
            if 'data_is_missing' not in dd.keys():
                pd[sid] = {'OAligned': 0, 'SAligned': 0, 'UAligned': 0}
                for mapq, cnt in dd.items():
                    if mapq == 'unmapped':
                        pd[sid]['UAligned'] += int(cnt)
                    elif int(mapq) >= 40:
                        pd[sid]['OAligned'] += int(cnt)
                    else:
                        pd[sid]['SAligned'] += int(cnt)

        self.add_section(
            name = 'Mapping Summary',
            anchor = 'biscuit-mapping',
            description = 'Shows the number of optimally aligned reads (defined by MAPQ >= 40), ' \
            'suboptimally aligned reads (MAPQ < 40), and unmapped reads. ' \
            'A good library should have a high fraction of reads that are optimally aligned. ' \
            'Note, Suboptimally aligned reads include both nonunique alignments and imperfect alignments.',
            plot = bargraph.plot(pd,
                                 OrderedDict([('OAligned', {'name': 'Optimally Aligned Reads'}),
                                              ('SAligned', {'name': 'Suboptimally Aligned Reads'}),
                                              ('UAligned', {'name': 'Unaligned Reads'})]),
                                 {'id': 'biscuit_mapping_summary',
                                  'title': 'BISCUIT: Mapping Summary',
                                  'ylab': 'Number of Reads',
                                  'cpswitch_counts_label': '# Reads'})
        )

        # Mapping quality distribution
        total = {}
        for sid, dd in self.mdata['align_mapq'].items():
            if 'data_is_missing' not in dd.keys():
                total[sid] = sum([int(cnt) for _, cnt in dd.items() if _ != 'unmapped'])

        pd_mapping = {}
        for sid, dd in self.mdata['align_mapq'].items():
            if 'data_is_missing' not in dd.keys():
                mapqcnts = []
                for mapq in range(61):
                    if str(mapq) in dd:
                        mapqcnts.append(100 * float(dd[str(mapq)]) / total[sid])
                    else:
                        mapqcnts.append(0)
                pd_mapping[sid] = dict(zip(range(61), mapqcnts))

        self.add_section(
            name = 'Mapping Quality Distribution',
            anchor = 'biscuit-mapq',
            description = 'Shows the distribution of primary mapping qualities as a percentage of the total number of mapped reads.',
            plot = linegraph.plot(pd_mapping,
                                  {'id': 'biscuit_mapping',
                                   'title': 'BISCUIT: Mapping Information',
                                   'ymin': 0,
                                   'yLabelFormat': '{value}%',
                                   'tt_label': '<strong>Q{point.x}:</strong> {point.y:.2f}% of reads',
                                   'name': 'Mapping Quality',
                                   'ylab': '% Primary Mapped Reads',
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
            log.warning('Missing data in {}. Will not fill corresponding entries.'.format(fn))
            return {'data_is_missing': 1}

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
            if 'data_is_missing' not in dd.keys():
                pd1[sid] = dd['read1']
                pd2[sid] = dd['read2']

        self.add_section(
            name = 'Mapping Strand Distribution',
            anchor = 'biscuit-strands',
            description = 'Shows the number of reads mapped to each strand. ' \
            'Most bisulfite libraries map read 1 to the parent strand (`++` or `--`) and ' \
            'read 2 to the daughter/synthesized strand (`+-` or `-+`). ' \
            'PBAT and most single-cell/low input libraries often do not follow this assumption.',
            plot = bargraph.plot([pd1, pd2],
                                 [OrderedDict([('++', {'name': '++: Waston-Aligned, Waston-Bisulfite Conversion', 'color': '#F53855'}),
                                               ('+-', {'name': '+-: Waston-Aligned, Crick-Bisulfite Conversion' , 'color': '#E37B40'}),
                                               ('-+', {'name': '-+: Crick-Aligned, Waston-Bisulfite Conversion' , 'color': '#46B29D'}),
                                               ('--', {'name': '--: Crick-Aligned, Crick-Bisulfite Conversion'  , 'color': '#324D5C'})]),
                                  OrderedDict([('++', {'name': '++: Waston-Aligned, Waston-Bisulfite Conversion', 'color': '#F53855'}),
                                               ('+-', {'name': '+-: Waston-Aligned, Crick-Bisulfite Conversion' , 'color': '#E37B40'}),
                                               ('-+', {'name': '-+: Crick-Aligned, Waston-Bisulfite Conversion' , 'color': '#46B29D'}),
                                               ('--', {'name': '--: Crick-Aligned, Crick-Bisulfite Conversion'  , 'color': '#324D5C'})])],
                                 {'id': 'biscuit_strands',
                                  'title': 'BISCUIT: Mapping Strand Distribution',
                                  'ylab': 'Number of Reads',
                                  'cpswitch_c_active': True,
                                  'cpswitch_counts_label': '# Reads',
                                  'data_labels': [{'name': 'Read 1', },
                                                  {'name': 'Read 2', }]})
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
            log.warning('Missing data in {}. Will not fill corresponding entries.'.format(fn))
            return {'data_is_missing': 1}

        data = {'I': {}, 'S': {}} # I = insert size, S = SW-score
        for l in file_data:
            fields = l.split()
            if fields[0] == 'I':
                data[fields[0]][int(fields[1])] = float(fields[2]) * 100
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
            if 'data_is_missing' not in dd.keys():
                pd_isize[sid] = dd['I']

        self.add_section(
            name = 'Insert Size Distribution',
            anchor = 'biscuit-isize',
            description = "Shows the distribution of insert sizes.",
            plot = linegraph.plot(pd_isize,
                                  {'id': 'biscuit_isize',
                                   'title': 'BISCUIT: Insert Size Distribution',
                                   'ymin': 0,
                                   'yLabelFormat': '{value}%',
                                   'smooth_points': 500, # limit number of points / smooth data
                                   'tt_label': '<strong>Q{point.x}:</strong> {point.y:.2f}% of reads',
                                   'ylab': '% Mapped Reads',
                                   'xlab': 'Insert Size'})
        )

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
            log.warning('Missing data in {}. Will not fill corresponding entries.'.format(fn))
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
            No returns, generates Cumulative Base Coverage chart
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
            name = 'Cumulative Base Coverage',
            anchor = 'biscuit-coverage-base',
            description = 'Shows the cumulative base and CpG coverage. ' \
            '"Q40" shows only those reads with mapping quality greater than or equal to 40. ' \
            'High and low GC content regions are the top and bottom 10% of 100bp windows for GC content, respectively.',
            plot = linegraph.plot(mdata,
                                  {'id': 'biscuit_coverage_base',
                                   'title': 'BISCUIT: Cumulative Base Coverage',
                                   'ymin': 0,
                                   'xLabelFormat': '{value}X',
                                   'xlab': 'Sequencing Depth',
                                   'data_labels': [{'name': 'All', 'ylab':'Millions of Bases'},
                                                   {'name': 'Q40', 'ylab':'Millions of Bases (Q40)'},
                                                   {'name': 'Q40 low GC', 'ylab':'Millions of Low-GC Bases (Q40)'},
                                                   {'name': 'Q40 high GC', 'ylab':'Millions of High-GC Bases (Q40)'},
                                                   {'name': 'CpG (all)', 'ylab':'Millions of CpGs'},
                                                   {'name': 'CpG Q40', 'ylab':'Millions of CpGs (Q40)'},
                                                   {'name': 'CpG Q40 low GC', 'ylab':'Millions of Low-GC CpGs (Q40)'},
                                                   {'name': 'CpG Q40 high GC', 'ylab':'Millions of High-GC CpGs (Q40)'}]})
        )

        basecov = OrderedDict()
        for sid, dd in mdata[0].items():
            if sid not in basecov:
                basecov[sid] = {}

            if dd[0] == -1:
                basecov[sid]['all'] = -1
            else:
                basecov[sid]['all'] = float(dd[1])/float(dd[0])*100.0

        for sid, dd in mdata[1].items():
            if sid not in basecov:
                basecov[sid] = {}

            if dd[0] == -1:
                basecov[sid]['uniq'] = -1
            else:
                basecov[sid]['uniq'] = float(dd[1])/float(dd[0])*100.0

        for sid, dd in mdata[4].items():
            if sid not in basecov:
                basecov[sid] = {}

            if dd[0] == -1:
                basecov[sid]['cpg_all'] = -1
            else:
                basecov[sid]['cpg_all'] = float(dd[1])/float(dd[0])*100.0

        for sid, dd in mdata[5].items():
            if sid not in basecov:
                basecov[sid] = {}

            if dd[0] == -1:
                basecov[sid]['cpg_uniq'] = -1
            else:
                basecov[sid]['cpg_uniq'] = float(dd[1])/float(dd[0])*100.0

        # base coverage >=1x table
        if len(basecov) == 0:
            return

        self.add_section(
            name = 'Coverage by At Least One Read',
            anchor = 'biscuit-coverage-base-table',
            description = 'The fraction of the genome and genomic CpGs covered by at least one read.',
            plot = table.plot(basecov,
                              {'all': {'title': 'Genome (All)', 'max': 100, 'min': 0, 'suffix': '%'},
                               'uniq': {'title': 'Genome (Unique)', 'max': 100, 'min': 0, 'suffix': '%'},
                               'cpg_all': {'title': 'CpG (All)', 'max': 100, 'min': 0, 'suffix': '%'},
                               'cpg_uniq': {'title': 'CpG (Unique)', 'max': 100, 'min': 0, 'suffix': '%'}})
        )

    def chart_covdist_q40(self):
        pass # handled in chart_covdist

    def chart_covdist_q40_botgc(self):
        pass # handled in chart_covdist

    def chart_covdist_q40_topgc(self):
        pass # handled in chart_covdist

    def chart_covdist_cpg(self):
        pass # handled in chart_covdist

    def chart_covdist_cpg_q40(self):
        pass # handled in chart_covdist

    def chart_covdist_cpg_q40_botgc(self):
        pass # handled in chart_covdist

    def chart_covdist_cpg_q40_topgc(self):
        pass # handled in chart_covdist

    ########################################
    ####      Depths and Uniformity     ####
    ########################################
