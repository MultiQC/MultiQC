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
            'qc_cv': {},
            'qc_cpg_cv': {},
            # CpG distribution
            'qc_cpg_dist': {},
            # Duplicate reporting
            'dup_report': {},
            'markdup_report': {},
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
            '_markdup_report',
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

        for sid, dd in self.mdata['dup_report'].items():
            if sid not in pd:
                pd[sid] = {}
            if 'all' in dd and dd['all'] is not None:
                pd[sid]['%dupRate_All'] = dd['all']
            if 'all-q40' in dd and dd['all-q40'] is not None:
                pd[sid]['%dupRate_All-q40'] = dd['all-q40']

        pheader = OrderedDict()
        pheader['%aligned']         = {'title':'% Aligned', 'max':100, 'min':0, 'suffix':'%', 'scale':'Greens'}
        pheader['%dupRate_All']     = {'title':'% Overall Dup Rate', 'max':100, 'min':0, 'suffix':'%', 'scale':'Reds'}
        pheader['%dupRate_All-q40'] = {'title':'% Q40 Overall Dup Rate', 'max':100, 'min':0, 'suffix':'%', 'scale':'Reds'}
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
            description = 'Shows the number of optimally aligned reads (defined by MAPQ>=40), ' \
            'suboptimally aligned reads (MAPQ<40), and unmapped reads. ' \
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
                    pd_wg[sid]['mu_'+ctg] = dd[ctg]['mu']
                    pd_wg[sid]['sigma_'+ctg] = dd[ctg]['sigma']
                    pd_wg[sid]['cv_'+ctg] = dd[ctg]['cv']

        pheader_wg = OrderedDict()
        pheader_wg['mu_all'] = {'title': 'Mu.Genome', 'description': 'Whole Genome Mean'}
        pheader_wg['mu_all_topgc'] = {'title': 'Mu.High.GC', 'description': 'Top Decile for GC Content Mean'}
        pheader_wg['mu_all_botgc'] = {'title': 'Mu.Low.GC', 'description': 'Bottom Decile for GC Content Mean'}
        pheader_wg['sigma_all'] = {'title': 'Sigma.Genome', 'description': 'Whole Genome Std. Dev.'}
        pheader_wg['sigma_all_topgc'] = {'title': 'Sigma.High.GC', 'description': 'Top Decile for GC Content Std. Dev.'}
        pheader_wg['sigma_all_botgc'] = {'title': 'Sigma.Low.GC', 'description': 'Bottom Decile for GC Content Std. Dev.'}
        pheader_wg['cv_all'] = {'title': 'CV.Genome', 'description': 'Whole Genome Coeff. of Var.'}
        pheader_wg['cv_all_topgc'] = {'title': 'CV.High.GC', 'description': 'Top Decile for GC Content Coeff. of Var.'}
        pheader_wg['cv_all_botgc'] = {'title': 'CV.Low.GC', 'description': 'Bottom Decile for GC Content Coeff. of Var.'}

        self.add_section(
            name = 'Sequencing Depth - Whole Genome',
            anchor = 'biscuit-seq-depth-all',
            description = 'Shows the sequence depth mean, standard deviation, and uniformity measured by the ' \
            'Coefficient of Variation (sigma/mu) for reads with MAPQ>=40 across the entire genome. The top and ' \
            'bottom deciles for GC content were measured on 100bp non-overlapping windows.',
            plot = table.plot(pd_wg, pheader_wg)
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
                    pd_cg[sid]['mu_'+ctg] = dd[ctg]['mu']
                    pd_cg[sid]['sigma_'+ctg] = dd[ctg]['sigma']
                    pd_cg[sid]['cv_'+ctg] = dd[ctg]['cv']

        pheader_cg = OrderedDict()
        pheader_cg['mu_cpg'] = {'title': 'CG.Mu.Genome', 'description': 'All CpGs Mean'}
        pheader_cg['mu_cpg_topgc'] = {'title': 'CG.Mu.High.GC', 'description': 'CpGs in Top Decile for GC Content Mean'}
        pheader_cg['mu_cpg_botgc'] = {'title': 'CG.Mu.Low.GC', 'description': 'CpGs in Bottom Decile for GC Content Mean'}
        pheader_cg['sigma_cpg'] = {'title': 'CG.Sigma.Genome', 'description': 'All CpGs Std. Dev.'}
        pheader_cg['sigma_cpg_topgc'] = {'title': 'CG.Sigma.High.GC', 'description': 'CpGs in Top Decile for GC Content Std. Dev.'}
        pheader_cg['sigma_cpg_botgc'] = {'title': 'CG.Sigma.Low.GC', 'description': 'CpGs in Bottom Decile for GC Content Std. Dev.'}
        pheader_cg['cv_cpg'] = {'title': 'CG.CV.Genome', 'description': 'All CpGs Coeff. of Var.'}
        pheader_cg['cv_cpg_topgc'] = {'title': 'CG.CV.High.GC', 'description': 'CpGs in Top Decile for GC Content Coeff. of Var.'}
        pheader_cg['cv_cpg_botgc'] = {'title': 'CG.CV.Low.GC', 'description': 'CpGs in Bottom Decile for GC Content Coeff. of Var.'}

        self.add_section(
            name = 'Sequencing Depth - CpGs Only',
            anchor = 'biscuit-seq-depth-cpg',
            description = 'Shows the sequence depth mean, standard deviation, and uniformity measured by the ' \
            'Coefficient of Variation (sigma/mu) at CpG locations in reads with MAPQ>=40. The top and ' \
            'bottom deciles for GC content were measured on 100bp non-overlapping windows.',
            plot = table.plot(pd_cg, pheader_cg)
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

        #TODO: Assumes this value is in _cpg_dist_table.txt -
        #      may need to change this in future
        m = re.search(r'#CpG Islands\t(\d+)', f, re.MULTILINE)
        num_cgi = -999
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
                data['cgi_coverage'][k] = 100 * float(m.group(1)) / num_cgi

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
        #    pd['Genome'][ctg1] = 100 * float(dd[ctg]['uc']) / dd['TotalCpGs']['uc']

        for sid, dd in self.mdata['qc_cpg_dist'].items():
            pd[sid] = OrderedDict()
            for ctg in ['ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
                ctg1 = ctg.replace('CpGs','')
                pd[sid][ctg1] = 100 * float(dd[ctg]['uc']) / dd['TotalCpGs']['uc']

        self.add_section(
            name = 'CpG Coverage by Genomic Feature',
            anchor = 'biscuit-coverage-cpg-dist',
            description = 'Shows the fraction of uniquely covered CpGs for different categories ' \
            '(exonic, repeat-masked, genic, and CpG islands) relative to the total number of uniquely ' \
            'covered CpGs in the dataset. CpGs are from reads with MAPQ>=40.',
            plot = table.plot(pd,
                              OrderedDict([('Exonic', {'title': 'Exonic CpGs', 'max': 100, 'min': 0, 'suffix': '%'}),
                                           ('Repeat', {'title': 'Repeat-Masked CpGs', 'max': 100, 'min': 0, 'suffix': '%'}),
                                           ('Genic', {'title': 'Genic CpGs', 'max': 100, 'min': 0, 'suffix': '%'}),
                                           ('CGI', {'title': 'CpG Island CpGs', 'max': 100, 'min': 0, 'suffix': '%'})]),
                              {'id': 'cpg-cov-gen-feature'})
        )

        # CpG Islands
        pd = dict([(sid, dd['cgi_coverage']) for sid, dd in self.mdata['qc_cpg_dist'].items()])

        self.add_section(
            name = 'CpG Island Coverage',
            anchor = 'biscuit-coverage-cgi',
            description = 'Shows the percentage of CpG islands (out of all CpG islands in the genome) ' \
            'that have at least 1, 3, 5, or 10 different CpGs covered. Coverage is based on reads with MAPQ>=40.',
            plot = table.plot(pd,
                              OrderedDict([('one', {'title': '>=1', 'suffix': '%', 'description': 'CpG islands with at least one CpG covered'}),
                                           ('three', {'title': '>=3', 'suffix': '%', 'description': 'CpG islands with at least three CpGs covered'}),
                                           ('five', {'title': '>=5', 'suffix': '%', 'description': 'CpG islands with at least five CpGs covered'}),
                                           ('ten', {'title': '>=10', 'suffix': '%', 'description': 'CpG islands with at least ten CpGs covered'})]),
                              {'id': 'cgi-cov-table'})
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
                data[k] = (100 * float(m2.group(1)) / float(m1.group(1)))
            else:
                data[k] = -1.0

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

        pd = dict([(sid, dd) for sid, dd in self.mdata['dup_report'].items()])

        self.add_section(
            name = 'Duplicate Rates',
            anchor = 'biscuit-dup-report',
            description = 'Shows the percentage of bases that are duplicates out of the total number of bases' \
            'High and low GC content regions are the top and bottom 10% of 100bp windows for GC content, respectively.',
            plot = table.plot(pd,
                              OrderedDict([('all', {'title': 'Overall', 'suffix': '%', 'max': 100, 'min': 0}),
                                           ('topGC', {'title': 'Overall High GC', 'suffix': '%', 'max': 100, 'min': 0}),
                                           ('lowGC', {'title': 'Overall Low GC', 'suffix': '%', 'max': 100, 'min': 0}),
                                           ('all-q40', {'title': 'Q40', 'suffix': '%', 'max': 100, 'min': 0}),
                                           ('topGC-q40', {'title': 'Q40 High GC', 'suffix': '%', 'max': 100, 'min': 0}),
                                           ('botGC-q40', {'title': 'Q40 Low GC', 'suffix': '%', 'max': 100, 'min': 0})]),
                              {'id': 'biscuit_dup_report'})
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
                data['dupRatePE'] = 100 * float(data['PEdup']) / data['PE']
            else:
                data['dupRatePE'] = None

        m = re.search(r'marked (\d+) duplicates from (\d+) single-end reads', f, re.MULTILINE)
        if m is None:
            data = {'SE': None, 'SEdup': None, 'dupRateSE': None}
        else:
            data['SEdup'] = int(m.group(1))
            data['SE'] = int(m.group(2))
            if data['SE'] > 0:
                data['dupRateSE'] = 100 * float(data['SEdup']) / data['SE']
            else:
                data['dupRateSE'] = None
                
        m = re.search(r'identified (\d+) dangling paired-end reads', f, re.MULTILINE)
        if m is None:
            data['PEdangling'] = None
        else:
            data['PEdangling'] = int(m.group(1))
        
        return data

    def chart_markdup_report(self):
        '''
        biscuit markdup is being deprecated, so pass until officially removed
        '''
        pass

