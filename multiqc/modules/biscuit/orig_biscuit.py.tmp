#!/usr/bin/env python

""" MultiQC module to parse output from BISCUIT """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc.plots import linegraph, bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

# Log parsing regexes
class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='BISCUIT', anchor='biscuit',
            href="https://github.com/zwdzwd/biscuit",
            info="is a tool to map bisulfite converted sequence reads and determine"\
            " cytosine methylation states.")

        # Set up data structures
        self.mdata = {
            'align_mapq': {},
            'align_strand': {},
            'align_isize': {},
            # cpg coverage and base coverage
            'covdist': {},
            'covdist_q40': {},
            'covdist_q40_botgc': {},
            'covdist_q40_topgc': {},
            'covdist_cpg': {},
            'covdist_cpg_q40': {},
            'covdist_cpg_q40_botgc': {},
            'covdist_cpg_q40_topgc': {},
            # uniformity
            'qc_cv': {},
            'qc_cpg_cv': {},
            # cpg distribtuion
            'qc_cpg_dist': {},
            # dup report
            'dup_report': {},
            'markdup_report': {},
            # retention
            'retention_dist': {},
            'retention_dist_byread': {},
            'retention_cpg_readpos': {},
            'retention_cph_readpos': {},
            'retention_rate_byread': {},
            'retention_rate_bybase': {},
        }

        file_suffixes = [
            '_mapq_table',
            '_strand_table',
            '_markdup_report',
            '_markdup.bam',
            '_markdup',
            '_isize_score_table',
            '_covdist_table',
            '_covdist_q40_table',
            '_covdist_q40_botgc_table',
            '_covdist_q40_topgc_table',
            '_covdist_cpg_table',
            '_covdist_cpg_q40_table',
            '_covdist_cpg_q40_botgc_table',
            '_covdist_cpg_q40_topgc_table',
            '_all_cv_table',
            '_cpg_cv_table',
            '_cpg_dist_table',
            '_beta_table',
            '_freqOfTotalRetentionPerRead',
            '_CpHRetentionByReadPos',
            '_CpGRetentionByReadPos',
            '_totalReadConversionRate']

        # Find and parse alignment reports
        for k in self.mdata:
            for f in self.find_log_files('biscuit/{}'.format(k)):

                # this cleans s_name before further processing
                s_name = self.clean_s_name(f['s_name'], f['root'])

                # clean file uncommon suffixes
                for file_suffix in file_suffixes:
                    s_name = s_name.replace(file_suffix, '')

                if s_name in self.mdata[k]:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))

                self.mdata[k][s_name] = getattr(self, 'parse_logs_%s' % k)(f['f'], f['fn'])

        for k in self.mdata:
            self.mdata[k] = self.ignore_samples(self.mdata[k])

        if sum([len(self.mdata[k]) for k in self.mdata]) == 0:
            raise UserWarning

        # Basic Stats Table
        self.biscuit_stats_table()

        # Write out to the report
        log.info("Found {} reports".format(sum([len(self.mdata[k]) for k in self.mdata])))
        for k in self.mdata:
            if len(self.mdata[k]) > 0:
                self.write_data_file(self.mdata[k], 'multiqc_biscuit_{}'.format(k))
                log.debug("Found %d %s reports" % (len(self.mdata[k]), k))
                getattr(self, 'chart_{}'.format(k))()

    def biscuit_stats_table(self):

        pd = {}
        for sid, dd in self.mdata['align_mapq'].items():
            allreads = sum([int(_) for _ in dd.values()])
            pd[sid] = {'%aligned': float(allreads-int(dd['unmapped']))/allreads*100}
        
        for sid, dd in self.mdata['markdup_report'].items():
            if sid not in pd:
                pd[sid] = {}
            if 'dupRatePE' in dd and dd['dupRatePE'] is not None:
                pd[sid]['%dupratePE'] = dd['dupRatePE']
            if 'dupRateSE' in dd and dd['dupRateSE'] is not None:
                pd[sid]['%duprateSE'] = dd['dupRateSE']

        pheader = OrderedDict()
        pheader['%aligned'] = {'title':'% Aligned', 'max':100, 'min':0, 'suffix':'%','scale':'Greens'}
        pheader['%dupratePE'] = {'title':'% DupsPE', 'max':100, 'min':0, 'suffix':'%','scale':'Reds'}
        pheader['%duprateSE'] = {'title':'% DupsSE', 'max':100, 'min':0, 'suffix':'%','scale':'Reds'}
        self.general_stats_addcols(pd, pheader)

    def parse_logs_align_mapq(self, f, fn): # _mapq_table.txt

        data = {}
        for l in f.splitlines()[2:]:
            s = l.split()
            data[s[0]] = s[1]   # mapping quality > number of reads

        return data

    def chart_align_mapq(self):

        # fraction of optimally mapped reads
        pd = {}
        for sid, dd in self.mdata['align_mapq'].items():
            pd[sid] = {'OAligned':0, 'SAligned':0, 'UAligned':1}
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
            description = 'This shows the fraction of optimally aligned reads, which is defined by mapQ >= 40.',
            helptext = 'A good library should have high fraction of reads optimally aligned. Suboptimally aligned reads include both nonunique alignments and imperfect alignments.',
            plot = bargraph.plot(pd, OrderedDict([
                ('OAligned', {'name':'Optimally Aligned Reads'}),
                ('SAligned', {'name':'Suboptimally Aligned Reads'}),
                ('UAligned', {'name':'Unaligned Reads'})
            ]), {'id':'biscuit_mapping_summary',
                 'title':'BISCUIT: Mapping Summary',
                 'ylab':'Number of Reads',
                 'cpswitch_counts_label': '# Reads'
            })
        )

        # Mapping quality together in one plot
        total = {}
        for sid, dd in self.mdata['align_mapq'].items():
            total[sid] = sum([int(cnt) for _, cnt in dd.items() if _ != "unmapped"])

        pd_mapping = {}
        for sid, dd in self.mdata['align_mapq'].items():
            mapqcnts = []
            for mapq in range(61):
                if str(mapq) in dd:
                    mapqcnts.append(float(dd[str(mapq)])/total[sid]*100)
                else:
                    mapqcnts.append(0)
            pd_mapping[sid] = dict(zip(range(61), mapqcnts))

        self.add_section(
            name = 'Mapping Quality Distribution',
            anchor = 'biscuit-mapq',
            description = "This plot shows the distribution of primary mapping quality.",
            plot = linegraph.plot(pd_mapping,
                {'id':'biscuit_mapping',
                 'title': 'BISCUIT: Mapping Information', 
                 'ymin': 0, 'yLabelFormat': '{value}%', 
                 'tt_label': '<strong>Q{point.x}:</strong> {point.y:.2f}% of reads',
                 'name':'Mapping Quality', 'ylab': '% Primary Mapped Reads','xlab': 'Mapping Quality'}))

    def parse_logs_align_strand(self, f, fn): # _strand_table.txt

        data = {'read1':{},'read2':{}}
        for l in f.splitlines()[1:]:
            m = re.search(r'strand\t([12])([+-]*)\t(\d+)', l)
            if m is not None:
                if m.group(1) == '1':
                    data['read1'][m.group(2)] = int(m.group(3))
                if m.group(1) == '2':
                    data['read2'][m.group(2)] = int(m.group(3))

        return data

    def chart_align_strand(self):

        # mapping strand distribution
        pd1 = {}
        pd2 = {}
        for sid, dd in self.mdata['align_strand'].items():
            pd1[sid] = dd['read1']
            pd2[sid] = dd['read2']

        self.add_section(
            name='Mapping Strand Distribution',
            anchor='biscuit-strands',
            description = "This plot shows the distribution of strand of mapping and strand of bisulfite conversion.",
            helptext="Most bisulfite libraries has read 1 goes to parent `++` or `--` and read 2 goes to daughter/synthesized `+-` or `-+`. PBAT or most single-cell/low input libraries typically don't observe this rule.",
            plot = bargraph.plot([pd1, pd2], 
                [OrderedDict([
                    ('++', {'name':'++: Waston-Aligned, Waston-Bisulfite Conversion', 'color': '#F53855'}),
                    ('+-', {'name':'+-: Waston-Aligned, Crick-Bisulfite Conversion', 'color': '#E37B40'}),
                    ('-+', {'name':'-+: Crick-Aligned, Waston-Bisulfite Conversion', 'color': '#46B29D'}),
                    ('--', {'name':'--: Crick-Aligned, Crick-Bisulfite Conversion', 'color': '#324D5C'}),]),
                OrderedDict([
                    ('++', {'name':'++: Waston-Aligned, Waston-Bisulfite Conversion', 'color': '#F53855'}),
                    ('+-', {'name':'+-: Waston-Aligned, Crick-Bisulfite Conversion', 'color': '#E37B40'}),
                    ('-+', {'name':'-+: Crick-Aligned, Waston-Bisulfite Conversion', 'color': '#46B29D'}),
                    ('--', {'name':'--: Crick-Aligned, Crick-Bisulfite Conversion', 'color': '#324D5C'})])], 
                {'id':'biscuit_strands',
                 'title':'BISCUIT: Mapping Strand Distribution',
                 'ylab':'Number of Reads',
                 'cpswitch_c_active': True,
                 'cpswitch_counts_label': '# Reads',
                 'data_labels': [
                    {'name': 'Read 1', },
                    {'name': 'Read 2', }]
            })
        )

    def parse_logs_align_isize(self, f, fn): # _isize_score_table.txt

        data = {}
        data['I'] = {}      # insert size
        data['S'] = {}      # SW-score

        for l in f.splitlines()[2:]:
            fields = l.split()
            if fields[0] == 'I':
                data[fields[0]][int(fields[1])] = float(fields[2]) * 100
            elif fields[0] == 'S':
                data[fields[0]][fields[1]] = float(fields[2])
                
        return data
    
    def chart_align_isize(self):

        pd_isize = {}
        for sid, dd in self.mdata['align_isize'].items():
            pd_isize[sid] = dd['I']

        self.add_section(
            name = 'Insert Size Distribution',
            anchor = 'biscuit-isize',
            description = "This plot shows the distribution of insert size.",
            plot = linegraph.plot(pd_isize, 
                {'id':'biscuit_isize', 'title': 'BISCUIT: Insert Size Distribution',
                 'ymin': 0, 'yLabelFormat': '{value}%', 
                 'smooth_points': 500,       # limit number of points / smooth data
                 'tt_label': '<strong>Q{point.x}:</strong> {point.y:.2f}% of reads', 
                 'ylab': '% Mapped Reads', 'xlab': 'Insert Size'}))

    #######################################
    #### CpG Coverage and Base Coverage ###
    #######################################
    def parse_logs_covdist(self, f, fn):

        data = {}
        ## handles the following tables:
        ## _covdist_q40_table.txt, _covdist_q40_botgc_table.txt
        ## _covdist_q40_topgc_table.txt, _covdist_table.txt
        ## _covdist_cpg_q40_table.txt, _covdist_cpg_q40_botgc_table.txt
        ## _covdist_cpg_q40_topgc_table.txt, _covdist_cpg_table.txt
        dd = {}
        for l in f.splitlines()[2:]:
            fields = l.split()
            dd[int(float(fields[0]))] = int(float(fields[1]))

        covs = sorted([k for k in dd])[:30]
        ccov_cnts = []
        _ccov_cnt = sum(dd.values())

        for cov in covs:
            ccov_cnts.append(_ccov_cnt/1000000.)
            _ccov_cnt -= dd[cov]
        data = dict(zip(covs,ccov_cnts))

        return data

    def parse_logs_covdist_q40(self, f, fn):
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_q40_botgc(self, f, fn):
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_q40_topgc(self, f, fn):
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_cpg(self, f, fn):
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_cpg_q40(self, f, fn):
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_cpg_q40_botgc(self, f, fn):
        return self.parse_logs_covdist(f, fn)

    def parse_logs_covdist_cpg_q40_topgc(self, f, fn):
        return self.parse_logs_covdist(f, fn)

    def chart_covdist(self):

        # base coverage
        basecov = OrderedDict()
        mdata = [
            self.mdata['covdist'],
            self.mdata['covdist_q40'],
            self.mdata['covdist_q40_botgc'],
            self.mdata['covdist_q40_topgc'],
            self.mdata['covdist_cpg'],
            self.mdata['covdist_cpg_q40'],
            self.mdata['covdist_cpg_q40_botgc'],
            self.mdata['covdist_cpg_q40_topgc']]

        self.add_section(
            name = 'Cumulative Base Coverage',
            anchor = 'biscuit-coverage-base',
            description = "This plot shows the cummulative base coverage. High and low GC content region are the top and bottom 10% 100bp window in GC content.",
            helptext = "Q40 means only reads mapped with mapping quality (Q) greater than or equal to 40 are considered.",
            plot = linegraph.plot(mdata, {'id':'biscuit_coverage_base',
                'title': 'BISCUIT: Cumulative Base Coverage',
                'xLabelFormat':'{value}X',
                'xlab': 'Sequencing Depth',
                'data_labels': [
                    {'name': 'All', 'ylab':'Million Bases'}, 
                    {'name': 'Q40', 'ylab':'Million Bases (Q40)'}, 
                    {'name': 'Q40 low GC', 'ylab':'Million Low-GC Bases (Q40)'}, 
                    {'name': 'Q40 high GC', 'ylab':'Million High-GC Bases (Q40)'}, 
                    {'name': 'CpG (all)', 'ylab':'Million CpGs'}, 
                    {'name': 'CpG Q40', 'ylab':'Million CpGs (Q40)'}, 
                    {'name': 'CpG Q40 low GC', 'ylab':'Million Low-GC CpGs (Q40)'}, 
                    {'name': 'CpG Q40 high GC', 'ylab':'Million High-GC CpGs (Q40)'}, 
                ]})
        )

        for sid, dd in mdata[0].items():
            if sid not in basecov:
                basecov[sid] = {}
            basecov[sid]['all'] = float(dd[1])/float(dd[0])*100.0

        for sid, dd in mdata[1].items():
            if sid not in basecov:
                basecov[sid] = {}
            basecov[sid]['uniq'] = float(dd[1])/float(dd[0])*100.0

        for sid, dd in mdata[4].items():
            if sid not in basecov:
                basecov[sid] = {}
            basecov[sid]['cpg_all'] = float(dd[1])/float(dd[0])*100.0

        for sid, dd in mdata[5].items():
            if sid not in basecov:
                basecov[sid] = {}
            basecov[sid]['cpg_uniq'] = float(dd[1])/float(dd[0])*100.0

        # base coverage >=1x table
        if len(basecov) == 0:
            return

        self.add_section(
            name = 'Coverage by At Least One Read',
            anchor = 'biscuit-coverage-base-table',
            description = 'The fraction of genome/genomic CpGs covered by at least one read.',
            plot = table.plot(basecov, {
                'all':{'title':'Genome (All)','max':100,'min':0,'suffix':'%'},
                'uniq':{'title':'Genome (Unique)','max':100,'min':0,'suffix':'%'},
                'cpg_all':{'title':'CpG (All)','max':100,'min':0,'suffix':'%'},
                'cpg_uniq':{'title':'CpG (Unique)','max':100,'min':0,'suffix':'%'},
            })
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

    ##############################
    #### Depths and Uniformity ###
    ##############################
    def parse_logs_qc_cv(self, f, fn): # _all_cv_table.txt

        data = {}
        targets = ['all','all_topgc','all_botgc']
        for k in targets:
            m = re.search('_{}\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)'.format(k), f, re.MULTILINE)
            if m is None:
                data[k] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data[k] = {
                    'mu': float(m.group(1)), 
                    'sigma': float(m.group(2)), 
                    'cv': float(m.group(3))}
        return data

    def chart_qc_cv(self): 

        # sequencing depth and uniformity
        pd = OrderedDict()
        for sid, dd in self.mdata['qc_cv'].items():
            pd[sid] = OrderedDict()
            for ctg in ['all','all_topgc','all_botgc']:
                if ctg in dd:
                    pd[sid]['cv_'+ctg] = dd[ctg]['cv']
                    pd[sid]['mu_'+ctg] = dd[ctg]['mu']

        for sid, dd in self.mdata['qc_cpg_cv'].items():
            if sid not in pd:
                pd[sid] = OrderedDict()
            for ctg in ['cpg','cpg_topgc','cpg_botgc']:
                if ctg in dd:
                    pd[sid]['cv_'+ctg] = dd[ctg]['cv']
                    pd[sid]['mu_'+ctg] = dd[ctg]['mu']

        pheader = OrderedDict()
        pheader['mu_all'] = {'title':'Mu.gnm','description':'Whole Genome'}
        pheader['mu_all_topgc'] = {'title':'Mu.high.gc','description':'Top Decile in GC Content'}
        pheader['mu_all_botgc'] = {'title':'Mu.low.gc','description':'Bottom Decile in GC Content'}
        pheader['cv_all'] = {'title':'CV.gnm','description':'Whole Genome'}
        pheader['cv_all_topgc'] = {'title':'CV.high.gc','description':'Top Decile in GC Content'}
        pheader['cv_all_botgc'] = {'title':'CV.low.gc','description':'Bottom Decile in GC Content'}
        pheader['mu_cpg'] = {'title':'CG.Mu.gnm','description':'All CpGs'}
        pheader['mu_cpg_topgc'] = {'title':'CG.Mu.high.gc','description':'Top Decile in GC Content'}
        pheader['mu_cpg_botgc'] = {'title':'CG.Mu.low.gc','description':'Bottom Decile in GC Content'}
        pheader['cv_cpg'] = {'title':'CG.CV.gnm','description':'All CpGs'}
        pheader['cv_cpg_topgc'] = {'title':'CG.CV.high.gc','description':'Top Decile in GC Content'}
        pheader['cv_cpg_botgc'] = {'title':'CG.CV.low.gc','description':'Bottom Decile in GC Content'}

        self.add_section(
            name = 'Sequencing Depth',
            anchor = 'biscuit-seq-depth',
            description = "This plot shows sequence depth mean and uniformity measured in Coefficient of Variation (mu/sigma), mapQ>40 only. CG.* shows the depth on CpG only. GC contents were measured on 100bp non-overlapping windows.",
            plot = table.plot(pd, pheader))

    def parse_logs_qc_cpg_cv(self, f, fn): # _cpg_cv_table.txt

        data = {}
        targets = ['cpg', 'cpg_topgc', 'cpg_botgc']
        for k in targets:
            m = re.search('_{}\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)'.format(k), f, re.MULTILINE)
            if m is None:
                data[k] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data[k] = {
                    'mu': float(m.group(1)), 
                    'sigma': float(m.group(2)), 
                    'cv': float(m.group(3))}
        return data

    def chart_qc_cpg_cv(self):
        pass # handled in chart_qc_cv

    #########################
    #### CpG Distribution ###
    #########################
    def parse_logs_qc_cpg_dist(self, f, fn): # _cpg_dist_table.txt

        data = {}
        for ctg in ['TotalCpGs', 'ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
            m = re.search('{}\t(\d+)\t(\d+)\t(\d+)'.format(ctg), f, re.MULTILINE)
            data[ctg] = {
                'a':int(m.group(1)),
                'uc':int(m.group(2)),
                'ac':int(m.group(3))}

        m = re.search(r'#CpG Islands\t(\d+)', f, re.MULTILINE)
        num_cgi = int(m.group(1))

        patterns = [
            (r'one CpG\t(\d+)', 'one'),
            (r'three CpGs\t(\d+)', 'three'),
            (r'five CpGs\t(\d+)', 'five'),
            (r'ten CpGs\t(\d+)', 'ten')]

        data['cgi_coverage'] = {}
        for pat, k in patterns:
            m = re.search(pat, f, re.MULTILINE)
            data['cgi_coverage'][k] = float(m.group(1))/num_cgi*100

        return data

    def chart_qc_cpg_dist(self): 

        # cpg distribution
        if len(self.mdata['qc_cpg_dist']) == 0:
            return 

        hdr = OrderedDict()
        pd = OrderedDict()
        sid = 'Genome'
        dd = list(self.mdata['qc_cpg_dist'].values())[0]
        pd[sid] = OrderedDict()
        for ctg in ['ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
            ctg1 = ctg.replace('CpGs','')
            hdr[ctg1] = {'max':100,'min':0,'suffix':'%'}
            pd[sid][ctg1] = float(dd[ctg]['uc']) / dd['TotalCpGs']['uc'] * 100

        hdr['Exonic']['description'] = 'Exonic CpGs'
        hdr['Repeat']['description'] = 'Repeat-Masked CpGs'
        hdr['Genic']['description']  = 'Genic CpGs'
        hdr['CGI']['description']    = 'CpG Island CpGs'
        for sid, dd in self.mdata['qc_cpg_dist'].items():
            pd[sid] = OrderedDict()
            for ctg in ['ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
                ctg1 = ctg.replace('CpGs','')
                pd[sid][ctg1] = float(dd[ctg]['uc']) / dd['TotalCpGs']['uc'] * 100

        self.add_section(
            name = 'CpG Coverage vs Genomic Features',
            anchor = 'biscuit-coverage-cpg-dist',
            description = "The top row shows how CpGs breaks down to different categories. Each other row shows the how CpGs uniquely covered by the given data breaks down to these categories. It is the fraction of CpGs in the given category out of all CpGs covered by the data.",
            plot = table.plot(pd, hdr)
        )

        pd = dict([(sid, dd['cgi_coverage']) for sid, dd in self.mdata['qc_cpg_dist'].items()])
        self.add_section(
            name = 'CpG Island Coverage',
            anchor = 'biscuit-coverage-cgi',
            description = "Each row shows the percentage of CpG islands (out of all CpG islands in the genome) that are covered in different numbers of CpGs. Coverage is based on reads with mapQ >= 40.",
            plot = table.plot(pd, OrderedDict([
                ('one', {'title':'>=1', 
                 'suffix':'%', 'description':'CpG islands with at least one CpG covered'}),
                ('three', {'title':'>=2', 
                 'suffix':'%', 'description':'CpG islands with at least three CpGs covered'}),
                ('five', {'title':'>=5', 
                 'suffix':'%', 'description':'CpG islands with at least five CpGs covered'}),
                ('ten', {'title':'>=10', 
                 'suffix':'%', 'description':'CpG islands with at least ten CpGs covered'}),
            ]), {'id':'cgi-cov-table'})
        )

    ###################
    #### dup report ###
    ###################

    def parse_logs_dup_report(self, f, fn): # _dup_report.txt

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
                data[k] = (float(m2.group(1)) / float(m1.group(1)) * 100)
            else:
                data[k] = 0.0

        return data

    def chart_dup_report(self):
        pass # duplication uniformity is currently skipped for brevity of the report

    def parse_logs_markdup_report(self, f, fn): # _markdup_report.txt

        data = {}
        m = re.search(r'marked (\d+) duplicates from (\d+) paired-end reads', f, re.MULTILINE)
        if m is None:
            data = {'PE':None,'PEdup':None,'dupRatePE':None}
        else:
            data['PE'] = int(m.group(2))
            data['PEdup'] = int(m.group(1))
            if data['PE'] > 0:
                data['dupRatePE'] = float(data['PEdup']) / data['PE'] * 100
            else:
                data['dupRatePE'] = None

        m = re.search(r'marked (\d+) duplicates from (\d+) single-end reads', f, re.MULTILINE)
        if m is None:
            data = {'SE':None,'SEdup':None,'dupRateSE':None}
        else:
            data['SE'] = int(m.group(2))
            data['SEdup'] = int(m.group(1))
            if data['SE'] > 0:
                data['dupRateSE'] = float(data['SEdup']) / data['SE'] * 100
            else:
                data['dupRateSE'] = None
                
        m = re.search(r'identified (\d+) dangling paired-end reads', f, re.MULTILINE)
        if m is None:
            data['PEdangling'] = None
        else:
            data['PEdangling'] = int(m.group(1))
        
        return data

    def chart_markdup_report(self):
        pass # duplication uniformity is currently skipped 

    ###################
    #### retention ####
    ###################

    def parse_logs_retention_dist(self, f, fn): # _CpGRetentionDist.txt
        
        sumcounts = 0
        for l in f.splitlines()[2:]:
            fields = l.split('\t')
            sumcounts += int(fields[1])
        
        data = {}
        for l in f.splitlines()[2:]:
            fields = l.split('\t')
            data[int(fields[0])] = int(fields[1]) / float(sumcounts)

        return data

    def chart_retention_dist(self):

        ## cytosine retention distribution
        mdata_meth = self.mdata['retention_dist']
        mdata = self.mdata['retention_dist_byread']

        pd = [
            mdata_meth,
            dict([(sid, dd['CA']) for sid, dd in mdata.items()]),
            dict([(sid, dd['CC']) for sid, dd in mdata.items()]),
            dict([(sid, dd['CG']) for sid, dd in mdata.items()]),
            dict([(sid, dd['CT']) for sid, dd in mdata.items()]),
        ]
        self.add_section(
            name = 'Number of Retention Distribution',
            anchor = 'biscuit-retention-read',
            description = "This plot shows the distribution of the number of retained cytosine in each read, up to 10.",
            plot = linegraph.plot(pd, {
                'id': 'biscuit_retention_read_cpa', 
                'xlab': 'Number of Retention within Read',
                'title': 'BISCUIT: Retention Distribution',
                'data_labels': [
                    {'name': 'CpG retention', 'ylab': 'Fraction of cytosine in CpG context', 'xlab': 'Retention Level (%)'},
                    {'name': 'Within-read CpA', 'ylab': 'Number of Reads'},
                    {'name': 'Within-read CpC', 'ylab': 'Number of Reads'},
                    {'name': 'Within-read CpG', 'ylab': 'Number of Reads'},
                    {'name': 'Within-read CpT', 'ylab': 'Number of Reads'},
                ]})
            )

    def parse_logs_retention_dist_byread(self, f, fn): # _freqOfTotalRetentionPerRead.txt

        data = OrderedDict([('CA',{}),('CC',{}),('CG',{}),('CT',{})])
        for l in f.splitlines()[2:]:
            fields = l.strip().split('\t')
            if fields[0] not in data:
                return {}
            if int(fields[1]) <= 10: # remove counts greater than 10
                data[fields[0]][int(fields[1])] = int(fields[2])

        return data

    def chart_retention_dist_byread(self):
        pass # handled in chart_retention_dist

    def parse_logs_retention_cpg_readpos(self, f, fn): # _CpGRetentionByReadPos.txt

        data = {}
        r1 = {'C':{}, 'R':{}}
        r2 = {'C':{}, 'R':{}}
        for l in f.splitlines()[2:]:
            fields = l.strip().split('\t')
            if fields[2] not in ['C','R'] or fields[0] not in ['1','2']:
                return {}
            if fields[0] == '1':
                r1[fields[2]][int(fields[1])] = int(fields[3])
            if fields[0] == '2':
                r2[fields[2]][int(fields[1])] = int(fields[3])

        r1rate = OrderedDict()
        for k in sorted(r1['C'].keys()):
            if k in r1['R']:
                r1rate[k] = float(r1['R'][k])/(r1['R'][k]+r1['C'][k])*100.0

        r2rate = {}
        for k in sorted(r2['C'].keys()):
            if k in r2['R']:
                r2rate[k] = float(r2['R'][k])/(r2['R'][k]+r2['C'][k])*100.0
            data = {'1':r1rate, '2':r2rate}

        return data

    def chart_retention_cpg_readpos(self):

        ## retention vs read position
        mdata = [
            dict([(k,v['1']) for k, v in self.mdata['retention_cph_readpos'].items()]),
            dict([(k,v['2']) for k, v in self.mdata['retention_cph_readpos'].items()]),
            dict([(k,v['1']) for k, v in self.mdata['retention_cpg_readpos'].items()]),
            dict([(k,v['2']) for k, v in self.mdata['retention_cpg_readpos'].items()]),
        ]
        self.add_section(
            name = 'Retention vs. Base Position in Read',
            anchor = 'biscuit-retention-cytosine',
            description = "This plot (aka. mbias plot) shows the distribution of cytosine retention rate in read.",
            plot = linegraph.plot(mdata, {
                'id': 'biscuit_retention_cytosine',
                'xlab': 'Position in Read', 'ymin':0, 'ymax':100, 'yMinRange':0, 'yFloor':0,
                'title': 'BISCUIT: Retention vs. Base Position in Read',
                'data_labels': [
                    {'name': 'CpH Read 1', 'ylab': 'CpH Retention Rate (%)', 'ymin':0, 'ymax':100},
                    {'name': 'CpH Read 2', 'ylab': 'CpH Retention Rate (%)', 'ymin':0, 'ymax':100},
                    {'name': 'CpG Read 1', 'ylab': 'CpG Retention Rate (%)', 'ymin':0, 'ymax':100},
                    {'name': 'CpG Read 2', 'ylab': 'CpG Retention Rate (%)', 'ymin':0, 'ymax':100},
                ]})
            )

    def parse_logs_retention_cph_readpos(self, f, fn):
        return self.parse_logs_retention_cpg_readpos(f, fn)

    def chart_retention_cph_readpos(self):
        pass # handled in chart_retention_cpg_readpos

    def parse_logs_retention_rate_byread(self, f, fn): # _totalReadConversionRate.txt

        data = {}
        m = re.match(r'([\d.]+)\t([\d.]+)\t([\d.]+)\t([\d.]+)', f.splitlines()[2])
        data['ca'] = float(m.group(1))*100.0
        data['cc'] = float(m.group(2))*100.0
        data['cg'] = float(m.group(3))*100.0
        data['ct'] = float(m.group(4))*100.0

        return data

    def chart_retention_rate_byread(self):

        mdata_byread = {}
        for sid, dd in self.mdata['retention_rate_byread'].items():
            sid = sid.replace('_totalReadConversionRate','')
            mdata_byread[sid] = dd

        mdata_bybase = {}
        for sid, dd in self.mdata['retention_rate_bybase'].items():
            sid = sid.replace('_totalBaseConversionRate','')
            mdata_bybase[sid] = dd

        pdata = {}
        for sid, dd in mdata_byread.items():
            pdata[sid] = dict(list(dd.items()) + list(mdata_bybase[sid].items()))

        pheader = OrderedDict()
        pheader['ca'] = {'title':'r.CpA', 'format':'{:,.2g}', 'min':0, 'max':100, 'description':'CpA', 'suffix':'%'}
        pheader['cc'] = {'title':'r.CpC', 'format':'{:,.2g}', 'min':0, 'max':100, 'description':'CpC', 'suffix':'%'}
        pheader['cg'] = {'title':'r.CpG', 'format':'{:,.2g}', 'min':0, 'max':100, 'description':'CpG', 'suffix':'%'}
        pheader['ct'] = {'title':'r.CpT', 'format':'{:,.2g}', 'min':0, 'max':100, 'description':'CpT', 'suffix':'%'}
        pheader['bca'] = {'title':'b.CpA', 'format':'{:,.2g}', 'min':0, 'max':100, 'description':'CpA', 'suffix':'%'}
        pheader['bcc'] = {'title':'b.CpC', 'format':'{:,.2g}', 'min':0, 'max':100, 'description':'CpC', 'suffix':'%'}
        pheader['bcg'] = {'title':'b.CpG', 'format':'{:,.2g}', 'min':0, 'max':100, 'description':'CpG', 'suffix':'%'}
        pheader['bct'] = {'title':'b.CpT', 'format':'{:,.2g}', 'min':0, 'max':100, 'description':'CpT', 'suffix':'%'}

        self.add_section(
            name = 'Cytosine Retention',
            anchor = 'biscuit-retention',
            description = "This plot shows cytosine retention rate. `r.` stands for read-averaging and `b.` stands for 'base-averaging.",
            helptext = "**Cytosine retention rate** is `1.0 - cytosine conversion rate`. Assuming full (complete but not over) bisulfite conversion, **cytosine retention rate** is the average cytosine modification (including 5mC, 5hmC etc) rate.",
            plot = table.plot(pdata, pheader))

    def parse_logs_retention_rate_bybase(self, f, fn): # _totalBaseConversionRate.txt

        data = {}
        m = re.match(r'([\d.]+)\t([\d.]+)\t([\d.]+)\t([\d.]+)', f.splitlines()[2])
        data['bca'] = float(m.group(1))*100.0
        data['bcc'] = float(m.group(2))*100.0
        data['bcg'] = float(m.group(3))*100.0
        data['bct'] = float(m.group(4))*100.0

        return data

    def chart_retention_rate_bybase(self):
        pass # handled in chart_retention_rate_byread



