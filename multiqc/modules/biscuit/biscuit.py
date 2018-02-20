#!/usr/bin/env python

""" MultiQC module to parse output from BISCUIT """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import beeswarm, linegraph, bargraph, table
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
            'beta': {},
            'mapq': {},
            'markdup': {},
            'coverage': {},
            'retention': {},
        }

        # Find and parse bismark alignment reports
        for k in self.mdata:
            for f in self.find_log_files('biscuit/%s' % k):
                d = getattr(self, 'parse_logs_%s' % k)(f['f'], f['fn'])
                self.mdata[k][f['fn']] = d

        num_parsed = 0
        for k in self.mdata:
            num_parsed += len(self.mdata[k])
        if num_parsed == 0:
            raise UserWarning

        # Basic Stats Table
        self.biscuit_stats_table()

        # Write out to the report
        for k in self.mdata:
            if len(self.mdata[k]) > 0:
                self.write_data_file(self.mdata[k], 'multiqc_biscuit_%s' % k)
                log.info("Found %d biscuit %s reports" % (len(self.mdata[k]), k))
                getattr(self, 'chart_biscuit_%s' % k)()

    def biscuit_stats_table(self):

        pd = {}
        for sid, dd in self.mdata['mapq'].iteritems():
            if not sid.endswith('_mapq_table.txt'):
                continue
            sid = sid.replace('_mapq_table.txt','')
            allreads = sum([int(_) for _ in dd.itervalues()])
            pd[sid] = {'%aligned':float(allreads-int(dd['unmapped']))/allreads*100}
        self.general_stats_addcols(pd, {'%aligned':{'title':'% Aligned', 'max':100, 'min':0, 'suffix':'%','scale':'Greens'}})

    def parse_logs_beta(self, f, fn):

        data = {}
        for l in f.splitlines():
            fields = l.split('\t')
            data[int(fields[0])] = int(fields[1])
            
        return data

    def chart_biscuit_beta(self):

        self.add_section(
            name = 'Beta Value Distribution',
            anchor = 'biscuit-beta',
            plot = linegraph.plot(self.mdata['beta'], {'ylab':'Count of cytosine in CpG context','xlab':'Beta Value (%)'})
        )
        
    def parse_logs_mapq(self, f, fn):

        data = {}
        if fn.endswith('_mapq_table.txt'):
            for l in f.splitlines():
                s = l.split()
                data[s[0]] = s[1]   # mapping quality > number of reads
                
        if fn.endswith('_strand_table.txt'):
            for l in f.splitlines():
                m = re.search(r'strand\t([+-]*)\t(\d+)', l)
                if m is not None:
                    data[m.group(1)] = int(m.group(2))
                
        return data
    
    def chart_biscuit_mapq(self):

        # fraction of optimally mapped reads
        pd = {}
        for sid, dd in self.mdata['mapq'].iteritems():
            if not sid.endswith('_mapq_table.txt'):
                continue
            sid = sid.replace('_mapq_table.txt','');
            pd[sid] = {'OAligned':0, 'SAligned':0, 'UAligned':1}
            for mapq, cnt in dd.iteritems():
                if mapq == 'unmapped':
                    pd[sid]['UAligned'] += int(cnt)
                elif int(mapq) >= 40:
                    pd[sid]['OAligned'] += int(cnt)
                else:
                    pd[sid]['SAligned'] += int(cnt)

        self.add_section(
            name = 'Mapping Efficiency',
            anchor = 'biscuit-mapq',
            plot = bargraph.plot(pd, OrderedDict([
                ('OAligned', {'name':'Optimally Aligned Reads', 'color': '#8bbc21'}),
                ('SAligned', {'name':'Suboptimally Aligned Reads', 'color': '#f7a35c'}),
                ('UAligned', {'name':'Unaligned Reads', 'color': '#000000'})
            ]), {'cpswitch_c_active': True})
        )

        # mapping quality distribution
        total = {}
        for sid, dd in self.mdata['mapq'].iteritems():
            if not sid.endswith('_mapq_table.txt'):
                continue
            sid = sid.replace('_mapq_table.txt','')
            total[sid] = sum([int(cnt) for _, cnt in dd.items() if _ != "unmapped"])

        pd = {}
        for sid, dd in self.mdata['mapq'].iteritems():
            if not sid.endswith('_mapq_table.txt'):
                continue
            sid = sid.replace('_mapq_table.txt','')
            mapqcnts = []
            for mapq in range(61):
                if str(mapq) in dd:
                    mapqcnts.append(float(dd[str(mapq)])/total[sid]*100)
                else:
                    mapqcnts.append(0)
            pd[sid] = dict(zip(range(61), mapqcnts))

        self.add_section(
            name = 'Mapping Quality',
            anchor = 'biscuit-mapq',
            description = "<p>This plot shows the distribution of Primary Mapping Quality.</p>",
            plot = linegraph.plot(pd, {'ylab': '% Primary Mapped Reads','xlab': 'Mapping Quality'})
        )

        # mapping strand distribution
        pd = dict([(k.replace('_strand_table.txt', ''),v) for k,v in self.mdata['mapq'].iteritems() if k.endswith('_strand_table.txt')])

        self.add_section(
            name='Mapping Strand Distribution',
            anchor='biscuit-mapq',
            description = "<p>This plot shows the distribution of strand of mapping and strand of bisulfite conversion.</p>",
            plot = bargraph.plot(pd, OrderedDict([
                ('++', {'name':'Waston-Aligned, Waston-Bisulfite Conversion', 'color': '#F53855'}),
                ('+-', {'name':'Waston-Aligned, Crick-Bisulfite Conversion', 'color': '#E37B40'}),
                ('-+', {'name':'Crick-Aligned, Waston-Bisulfite Conversion', 'color': '#46B29D'}),
                ('--', {'name':'Crick-Aligned, Crick-Bisulfite Conversion', 'color': '#324D5C'})
            ]), {'cpswitch_c_active': True})
        )

    def parse_logs_markdup(self, f, fn):

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

    def chart_biscuit_markdup(self):

        mdata = dict([(sid.replace('_markdup_report.txt',''), {'Duplication Rate':dd['dupRatePE']}) for sid, dd in self.mdata['markdup'].iteritems()])
        if len(mdata)>0 and all([dd['Duplication Rate'] is not None for dd in mdata.itervalues()]):
            self.add_section(
                name = 'Read Duplication Rate PE',
                anchor = 'biscuit-markdup',
                description = "<p>This plot shows paired-end (PE) read duplication rate.</p>",
                plot = table.plot(mdata, {'Duplication Rate':{'suffix':'%','min':0,'max':100}})
            )

        mdata = dict([(sid.replace('_markdup_report.txt',''), {'Duplication Rate':dd['dupRateSE']}) for sid, dd in self.mdata['markdup'].iteritems()])
        if len(mdata)>0 and all([dd['Duplication Rate'] is not None for dd in mdata.itervalues()]):
            self.add_section(
                name = 'Read Duplication Rate SE',
                anchor = 'biscuit-markdup',
                description = "<p>This plot shows single-end (SE) read duplication rate.</p>",
                plot = table.plot(mdata, {'Duplication Rate':{'suffix':'%','min':0,'max':100}})
            )

    def parse_logs_coverage(self, f, fn):

        data = {}
        if fn.endswith('_cpg_dist_table.txt'):
            
            for ctg in ['TotalCpGs', 'ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
                m = re.search(r'%s\t(\d+)\t(\d+)\t(\d+)' % ctg, f, re.MULTILINE)
                if m is None:
                    return None
                data[ctg] = {'a':int(m.group(1)),
                             'uc':int(m.group(2)),
                             'ac':int(m.group(3))}

            m = re.search(r'#CpG Islands\t(\d+)', f, re.MULTILINE)
            num_cgi = int(m.group(1))

            data['cgi_coverage'] = {}
            m = re.search(r'one CpG\t(\d+)', f, re.MULTILINE)
            data['cgi_coverage']['one'] = float(m.group(1))/num_cgi*100
            m = re.search(r'three CpGs\t(\d+)', f, re.MULTILINE)
            data['cgi_coverage']['three'] = float(m.group(1))/num_cgi*100
            m = re.search(r'five CpGs\t(\d+)', f, re.MULTILINE)
            data['cgi_coverage']['five'] = float(m.group(1))/num_cgi*100
            m = re.search(r'ten CpGs\t(\d+)', f, re.MULTILINE)
            data['cgi_coverage']['ten'] = float(m.group(1))/num_cgi*100

        elif fn.endswith('_cpg_cv_table.txt'):

            m = re.search(r'_cpg\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)', f, re.MULTILINE)
            if m is None:
                data['cpg'] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data['cpg'] = {'mu': float(m.group(1)), 'sigma': float(m.group(2)), 'cv': float(m.group(3))}
                
            m = re.search(r'_cpg_topgc\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)', f, re.MULTILINE)
            if m is None:
                data['cpg_topgc'] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data['cpg_topgc'] = {'mu': float(m.group(1)), 'sigma': float(m.group(2)), 'cv': float(m.group(3))}
            
            m = re.search(r'_cpg_botgc\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)', f, re.MULTILINE)
            if m is None:
                data['cpg_topgc'] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data['cpg_botgc'] = {'mu': float(m.group(1)), 'sigma': float(m.group(2)), 'cv': float(m.group(3))}
                
        elif fn.endswith('_all_cv_table.txt'):

            m = re.search(r'_all\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)', f, re.MULTILINE)
            if m is None:
                data['all'] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data['all'] = {'mu': float(m.group(1)), 'sigma': float(m.group(2)), 'cv': float(m.group(3))}
                
            m = re.search(r'_all_topgc\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)', f, re.MULTILINE)
            if m is None:
                data['all_topgc'] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data['all_topgc'] = {'mu': float(m.group(1)), 'sigma': float(m.group(2)), 'cv': float(m.group(3))}

            m = re.search(r'_all_botgc\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)', f, re.MULTILINE)
            if m is None:
                data['all_botgc'] = {'mu': -1, 'sigma': -1, 'cv': -1}
            else:
                data['all_botgc'] = {'mu': float(m.group(1)), 'sigma': float(m.group(2)), 'cv': float(m.group(3))}

        elif fn.endswith('_dup_report.txt'):

            m1 = re.search(r'#bases covered by all reads: (\d+)', f, re.MULTILINE)
            m2 = re.search(r'#bases covered by duplicate reads: (\d+)', f, re.MULTILINE)
            data['all'] = (float(m2.group(1)) / float(m1.group(1)) * 100) if (m1 is not None and m2 is not None and float(m1.group(1))>0) else 0.0

            m1 = re.search(r'#high-GC bases covered by all reads: (\d+)', f, re.MULTILINE)
            m2 = re.search(r'#high-GC bases covered by duplicate reads: (\d+)', f, re.MULTILINE)
            data['topGC'] = (float(m2.group(1)) / float(m1.group(1)) * 100) if (m1 is not None and m2 is not None and float(m1.group(1))>0) else 0.0
            
            m1 = re.search(r'#low-GC bases covered by all reads: (\d+)', f, re.MULTILINE)
            m2 = re.search(r'#low-GC bases covered by duplicate reads: (\d+)', f, re.MULTILINE)
            data['botGC'] = (float(m2.group(1)) / float(m1.group(1)) * 100) if (m1 is not None and m2 is not None and float(m1.group(1))>0) else 0.0

            m1 = re.search(r'#bases covered by all q40-reads: (\d+)', f, re.MULTILINE)
            m2 = re.search(r'#bases covered by duplicate q40-reads: (\d+)', f, re.MULTILINE)
            data['all-q40'] = (float(m2.group(1)) / float(m1.group(1)) * 100) if (m1 is not None and m2 is not None and float(m1.group(1))>0) else 0.0

            m1 = re.search(r'#high-GC bases covered by all q40-reads: (\d+)', f, re.MULTILINE)
            m2 = re.search(r'#high-GC bases covered by duplicate q40-reads: (\d+)', f, re.MULTILINE)
            data['topGC-q40'] = (float(m2.group(1)) / float(m1.group(1)) * 100) if (m1 is not None and m2 is not None and float(m1.group(1))>0) else 0.0
                
            m1 = re.search(r'#low-GC bases covered by all q40-reads: (\d+)', f, re.MULTILINE)
            m2 = re.search(r'#low-GC bases covered by duplicate q40-reads: (\d+)', f, re.MULTILINE)
            data['botGC-q40'] = (float(m2.group(1)) / float(m1.group(1)) * 100) if (m1 is not None and m2 is not None and float(m1.group(1))>0) else 0.0
                
        elif re.search(r'_covdist.*_table.txt', fn):

            ## handles the following tables:
            ## _covdist_q40_table.txt, _covdist_q40_botgc_table.txt
            ## _covdist_q40_topgc_table.txt, _covdist_table.txt
            ## _covdist_cpg_q40_table.txt, _covdist_cpg_q40_botgc_table.txt
            ## _covdist_cpg_q40_topgc_table.txt, _covdist_cpg_table.txt
            dd = {}
            for l in f.splitlines():
                fields = l.split()
                dd[int(float(fields[0]))] = int(float(fields[1]))

            covs = sorted([k for k in dd])[:30]
            ccov_cnts = []
            _ccov_cnt = sum(dd.values())
            for cov in covs:
                ccov_cnts.append(_ccov_cnt/1000000.)
                _ccov_cnt -= dd[cov]

            data = dict(zip(covs,ccov_cnts))
        else:
            raise Exception("Unknown file received: %s" % fn)
                
        return data

    def chart_biscuit_coverage(self):

        # base coverage
        basecov = OrderedDict()
        mdata = dict([(k.replace('_covdist_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_covdist_table.txt')])
        if len(mdata)>0:
            self.add_section(
                name = 'Cumulative Base Coverage',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows the cummulative base coverage.</p>",
                plot = linegraph.plot(mdata, {'ylab':'Million Bases'})
            )
            
            for sid, dd in mdata.iteritems():
                if sid not in basecov:
                    basecov[sid] = {}
                basecov[sid]['all'] = dd[1]/dd[0]*100

        # base coverage uniq
        mdata = dict([(k.replace('_covdist_q40_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_covdist_q40_table.txt')])
        if len(mdata)>0:
            self.add_section(
                name = 'Cumulative Base Coverage Q40',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows the cummulative base coverage, mapQ>40 only.</p>",
                plot = linegraph.plot(mdata, {'ylab':'Million Bases'})
            )

            for sid, dd in mdata.iteritems():
                if sid not in basecov:
                    basecov[sid] = {}
                basecov[sid]['uniq'] = dd[1]/dd[0]*100

        # base coverage >=1x table
        if len(basecov)>0:
            self.add_section(
                name = 'Base Coverage',
                anchor = 'biscuit-coverage',
                description = '<p>The fraction of genome covered by at least one read.</p>',
                plot = table.plot(basecov, {
                    'all':{'title':'All Reads','max':100,'min':0,'suffix':'%'},
                    'uniq':{'title':'Uniquely Mapped Reads','max':100,'min':0,'suffix':'%'},
                })
            )

        # base coverage low GC content
        mdata = dict([(k.replace('_covdist_q40_botgc_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_covdist_q40_botgc_table.txt')])
        if len(mdata)>0:
            self.add_section(
                name = 'Cumulative Base Coverage Q40, low GC content',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows the cummulative base coverage of low GC content region (bottom 10%), mapQ>40 only.</p>",
                plot = linegraph.plot(mdata, {'ylab':'Million Bases'})
            )

        # base coverage high GC content
        mdata = dict([(k.replace('_covdist_q40_topgc_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_covdist_q40_topgc_table.txt')])
        if len(mdata)>0:
            self.add_section(
                name = 'Cumulative Base Coverage Q40, high GC content',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows the cummulative base coverage of high GC content region (top 10%), mapQ>40 only.</p>",
                plot = linegraph.plot(mdata, {'ylab':'Million Bases'})
            )
        
        # cpg coverage
        cpgcov = OrderedDict()
        mdata = dict([(k.replace('_covdist_cpg_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_covdist_cpg_table.txt')])
        if len(mdata)>0:
            self.add_section(
                name = 'Cumulative CpG Coverage',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows the cummulative CpG coverage.</p>",
                plot = linegraph.plot(mdata, {'ylab':'Million CpGs'})
            )

            for sid, dd in mdata.iteritems():
                if sid not in cpgcov:
                    cpgcov[sid] = {}
                cpgcov[sid]['all'] = dd[1]/dd[0]*100            

        # cpg coverage uniq
        mdata = dict([(k.replace('_covdist_cpg_q40_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_covdist_cpg_q40_table.txt')])
        if len(mdata)>0:
            self.add_section(
                name = 'Cumulative CpG Coverage Q40',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows the cummulative CpG coverage, mapQ>40 only.</p>",
                plot = linegraph.plot(mdata, {'ylab':'Million CpGs'})
            )

            for sid, dd in mdata.iteritems():
                if sid not in cpgcov:
                    cpgcov[sid] = {}
                cpgcov[sid]['uniq'] = dd[1]/dd[0]*100

        # cpg coverage >= 1x table
        if len(cpgcov)>0:
            self.add_section(
                name = 'CpG Coverage',
                anchor = 'biscuit-coverage',
                description = '<p>The fraction of CpGs covered by at least one read.</p>',
                plot = table.plot(cpgcov, {
                    'all':{'title':'All Reads','max':100,'min':0,'suffix':'%'},
                    'uniq':{'title':'Uniquely Mapped Reads','max':100,'min':0,'suffix':'%'},
                })
            )

        # cpg coverage low GC content
        mdata = dict([(k.replace('_covdist_cpg_q40_botgc_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_covdist_cpg_q40_botgc_table.txt')])
        if len(mdata)>0:
            self.add_section(
                name = 'Cumulative CpG Coverage Q40, Low GC Content Region',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows the cummulative CpG coverage of low GC content region (bottom 10%), mapQ>40 only.</p>",
                plot = linegraph.plot(mdata, {'ylab':'Million Bases'})
            )

        # cpg coverage high GC content
        mdata = dict([(k.replace('_covdist_cpg_q40_topgc_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_covdist_cpg_q40_topgc_table.txt')])
        if len(mdata)>0:
            self.add_section(
                name = 'Cumulative Base Coverage Q40, High GC Content Region',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows the cummulative CpG coverage of high GC content region (top 10%), mapQ>40 only.</p>",
                plot = linegraph.plot(mdata, {'ylab':'Million Bases'})
            )        

        # cpg distribution
        mdata = dict([(k.replace('_cpg_dist_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_cpg_dist_table.txt')])
        if len(mdata)>0:
            hdr = OrderedDict()
            pd = OrderedDict()
            sid = 'Genome'
            dd = mdata.values()[0]
            pd[sid] = OrderedDict()
            for ctg in ['ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
                ctg1 = ctg.replace('CpGs','')
                hdr[ctg1] = {'max':100,'min':0,'suffix':'%'}
                pd[sid][ctg1] = float(dd[ctg]['uc']) / dd['TotalCpGs']['uc'] * 100

            hdr['Exonic']['description'] = 'Exonic CpGs'
            hdr['Repeat']['description'] = 'Repeat-Masked CpGs'
            hdr['Genic']['description']  = 'Genic CpGs'
            hdr['CGI']['description']    = 'CpG Island CpGs'
            for sid, dd in mdata.iteritems():
                pd[sid] = OrderedDict()
                for ctg in ['ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
                    ctg1 = ctg.replace('CpGs','')
                    pd[sid][ctg1] = float(dd[ctg]['uc']) / dd['TotalCpGs']['uc'] * 100

            self.add_section(
                name = 'CpG Coverage Distribution',
                anchor = 'biscuit-coverage',
                description = "<p>The top row shows how CpGs breaks down to different categories. Each other row shows the how CpGs uniquely covered by the given data breaks down to these categories. It is the fraction of CpGs in the given category out of all CpGs covered by the data.</p>",
                plot = table.plot(pd, hdr)
            )

            pd = dict([(sid, dd['cgi_coverage']) for sid, dd in mdata.iteritems()])
            self.add_section(
                name = 'CpG Island Coverage',
                anchor = 'biscuit-coverage',
                description = "<p>Each row shows the percentage of CpG islands (out of all CpG islands in the genome) that are covered in different numbers of CpGs. Coverage is based on reads with mapQ >= 40.</p>",
                plot = table.plot(pd, OrderedDict([
                    ('one', {'title':'>=1', 'suffix':'%',
                             'description':'CpG islands with at least one CpG covered'}),
                    ('three', {'title':'>=2', 'suffix':'%',
                               'description':'CpG islands with at least three CpGs covered'}),
                    ('five', {'title':'>=5', 'suffix':'%',
                              'description':'CpG islands with at least five CpGs covered'}),
                    ('ten', {'title':'>=10', 'suffix':'%',
                             'description':'CpG islands with at least ten CpGs covered'}),
                ]))
            )
        
        # base uniformity
        mdata = dict([(k.replace('_all_cv_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_all_cv_table.txt')])
        if len(mdata)>0:
            pd = OrderedDict()
            for sid, dd in mdata.iteritems():
                pd[sid] = OrderedDict()
                for ctg in ['all','all_topgc','all_botgc']:
                    if ctg in dd:
                        pd[sid][ctg] = dd[ctg]['cv']
            self.add_section(
                name = 'Sequence Depth Uniformity',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows sequence depth uniformity measured in Coefficient of Variation (mu/sigma), mapQ>40 only. GC contents were measured on 100bp non-overlapping windows.</p>",
                plot = table.plot(pd, OrderedDict([
                    ('all', {'title':'Genome','description':'Whole Genome Average'}),
                    ('all_topgc', {'title':'Highest GC','description':'Top Decile in GC Content'}),
                    ('all_botgc', {'title':'Lowest GC','description':'Bottom Decile in GC Content'})
                ]))
            )

            pd = OrderedDict()
            for sid, dd in mdata.iteritems():
                pd[sid] = OrderedDict()
                for ctg in ['all','all_topgc','all_botgc']:
                    if ctg in dd:
                        pd[sid][ctg] = dd[ctg]['mu']
            self.add_section(
                name = 'Sequence Depth Mean',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows sequence depth mean, mapQ>40 only.</p>",
                plot = table.plot(pd, OrderedDict([
                    ('all', {'title':'Genome','description':'Whole Genome Average'}),
                    ('all_topgc', {'title':'Highest GC','description':'Top Decile in GC Content'}),
                    ('all_botgc', {'title':'Lowest GC','description':'Bottom Decile in GC Content'})
                ]))
            )

            pd = OrderedDict()
            for sid, dd in mdata.iteritems():
                pd[sid] = OrderedDict()
                for ctg in ['all','all_topgc','all_botgc']:
                    if ctg in dd:
                        pd[sid][ctg] = dd[ctg]['sigma']
            self.add_section(
                name = 'Sequence Depth Sigma',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows sequence depth SD, mapQ>40 only.</p>",
                plot = table.plot(pd, OrderedDict([
                    ('all', {'title':'Genome','description':'Whole Genome Average'}),
                    ('all_topgc', {'title':'Highest GC','description':'Top Decile in GC Content'}),
                    ('all_botgc', {'title':'Lowest GC','description':'Bottom Decile in GC Content'})
                ]))
            )

        # cpg uniformity
        mdata = dict([(k.replace('_cpg_cv_table.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_cpg_cv_table.txt')])
        if len(mdata)>0:
            pd = OrderedDict()
            for sid, dd in mdata.iteritems():
                pd[sid] = OrderedDict()
                for ctg in ['cpg','cpg_topgc','cpg_botgc']:
                    if ctg in dd:
                        pd[sid][ctg] = dd[ctg]['cv']
            self.add_section(
                name = 'CpG Sequence Depth Uniformity',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows CpG sequence depth uniformity measured in Coefficient of Variation (mu/sigma), mapQ>40 only. GC contents were measured on 100bp non-overlapping windows.</p>",
                plot = table.plot(pd, OrderedDict([
                    ('cpg', {'title':'All CpGs','description':'All CpG Average'}),
                    ('cpg_topgc', {'title':'Highest GC','description':'Top Decile in GC Content'}),
                    ('cpg_botgc', {'title':'Lowest GC','description':'Bottom Decile in GC Content'})
                ]))
            )

            pd = OrderedDict()
            for sid, dd in mdata.iteritems():
                pd[sid] = OrderedDict()
                for ctg in ['cpg','cpg_topgc','cpg_botgc']:
                    if ctg in dd:
                        pd[sid][ctg] = dd[ctg]['mu']
            self.add_section(
                name = 'CpG Sequence Depth Mean',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows CpG sequence depth mean, mapQ>40 only.</p>",
                plot = table.plot(pd, OrderedDict([
                    ('cpg', {'title':'All CpGs','description':'All CpG Average'}),
                    ('cpg_topgc', {'title':'Highest GC','description':'Top Decile in GC Content'}),
                    ('cpg_botgc', {'title':'Lowest GC','description':'Bottom Decile in GC Content'})
                ]))
            )

            pd = OrderedDict()
            for sid, dd in mdata.iteritems():
                pd[sid] = OrderedDict()
                for ctg in ['cpg','cpg_topgc','cpg_botgc']:
                    if ctg in dd:
                        pd[sid][ctg] = dd[ctg]['sigma']
            self.add_section(
                name = 'CpG Sequence Depth Sigma',
                anchor = 'biscuit-coverage',
                description = "<p>This plot shows CpG sequence depth SD, mapQ>40 only.</p>",
                plot = table.plot(pd, OrderedDict([
                    ('cpg', {'title':'All CpGs','description':'All CpG Average'}),
                    ('cpg_topgc', {'title':'Highest GC','description':'Top Decile in GC Content'}),
                    ('cpg_botgc', {'title':'Lowest GC','description':'Bottom Decile in GC Content'})
                ]))
            )

        # duplication uniformity
        mdata = dict([(k.replace('_dup_report.txt',''),v) for k, v in self.mdata['coverage'].iteritems() if k.endswith('_dup_report.txt')])
        if len(mdata)>0:
            self.add_section(
                name = 'Base Coverage by Read Duplication',
                anchor = 'biscuit-coverage',
                description = "<p>This table shows base coverage by duplicate reads as a percentage of coverage by all reads.</p>",
                plot = table.plot(mdata, OrderedDict([
                    ('all', {'title':'Genome', 'suffix':'%', 'description':'Fraction of coverage by duplicate reads'}),
                    ('topGC', {'title':'HighGC', 'suffix':'%', 'description':'Fraction of coverage in top decile in GC content by duplicate reads'}),
                    ('botGC', {'title':'LowGC', 'suffix':'%', 'description':'Fraction of coverage in bottom decile in GC content by duplicate reads'}),
                    ('all-q40', {'title':'Genome-uniq', 'suffix':'%', 'description':'Fraction of coverage by duplicate reads with mapq>=40'}),
                    ('topGC-q40', {'title':'HighGC-uniq', 'suffix':'%', 'description':'Fraction of coverage in top decile in GC content by duplicate reads with mapq>=40'}),
                    ('botGC-q40', {'title':'LowGC-uniq', 'suffix':'%', 'description':'Fraction of coverage in bottom decile in GC content by duplicate reads with mapq>=40'}),
                ]))
            )

    def parse_logs_retention(self, f, fn):

        data = {}
        if fn.endswith('_totalBaseConversionRate.txt') or fn.endswith('_totalReadConversionRate.txt'):
            m = re.match(r'([\d.]+)\t([\d.]+)\t([\d.]+)\t([\d.]+)', f.splitlines()[1])
            data['ca'] = float(m.group(1))*100
            data['cc'] = float(m.group(2))*100
            data['cg'] = float(m.group(3))*100
            data['ct'] = float(m.group(4))*100

        elif fn.endswith('_freqOfTotalRetentionPerRead.txt'):

            data = OrderedDict([('CA',{}),('CC',{}),('CG',{}),('CT',{})])
            for i,l in enumerate(f.splitlines()):
                if i==0:
                    continue
                fields = l.strip().split('\t')
                if fields[0] not in data:
                    return {}
                if int(fields[1]) <= 10: # remove counts greater than 10
                    data[fields[0]][int(fields[1])] = int(fields[2])

        elif fn.endswith('_CpHRetentionByReadPos.txt'):
            r1 = {'C':{}, 'R':{}}
            r2 = {'C':{}, 'R':{}}
            for l in f.splitlines():
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
                    r1rate[k] = float(r1['R'][k])/(r1['R'][k]+r1['C'][k])*100
            r2rate = {}
            for k in sorted(r2['C'].keys()):
                if k in r2['R']:
                    r2rate[k] = float(r2['R'][k])/(r2['R'][k]+r2['C'][k])*100

            data = {'1':r1rate, '2':r2rate}
            
        return data

    def chart_biscuit_retention(self):

        mdata = dict([(k.replace('_totalReadConversionRate.txt',''),v) for k, v in self.mdata['retention'].iteritems() if k.endswith('_totalReadConversionRate.txt')])
        if len(mdata) > 0:
            pd = dict([(sid,dd) for sid, dd in mdata.iteritems()])
            self.add_section(
                name = 'Cytosine Retention',
                anchor = 'biscuit-retention',
                description = "<p>This plot shows cytosine retention rate by averaging retention level from each cytosine mapped in each read.</p>",
                plot = table.plot(pd, OrderedDict([
                    ('ca', {'title':'CpA','description':'CpA dinucleotide context','suffix':'%'}),
                    ('cc', {'title':'CpC','description':'CpC dinucleotide context','suffix':'%'}),
                    ('cg', {'title':'CpG','description':'CpG dinucleotide context','suffix':'%'}),
                    ('ct', {'title':'CpT','description':'CpT dinucleotide context','suffix':'%'})
                ]))
            )

        mdata = dict([(k.replace('_totalBaseConversionRate.txt',''),v) for k, v in self.mdata['retention'].iteritems() if k.endswith('_totalBaseConversionRate.txt')])
        if len(mdata) > 0:
            pd = dict([(sid,dd) for sid, dd in mdata.iteritems()])
            self.add_section(
                name = 'Base Averaged Cytosine Retention',
                anchor = 'biscuit-retention',
                description = "<p>This plot shows cytosine retention rate by averaging retention level from all cytosine bases.</p>",
                plot = table.plot(pd, OrderedDict([
                    ('ca', {'title':'CpA','description':'CpA dinucleotide context','suffix':'%'}),
                    ('cc', {'title':'CpC','description':'CpC dinucleotide context','suffix':'%'}),
                    ('cg', {'title':'CpG','description':'CpG dinucleotide context','suffix':'%'}),
                    ('ct', {'title':'CpT','description':'CpT dinucleotide context','suffix':'%'})
                ]))
            )

        mdata = dict([(k.replace('_freqOfTotalRetentionPerRead.txt',''),v) for k, v in self.mdata['retention'].iteritems() if k.endswith('_freqOfTotalRetentionPerRead.txt')])
        if len(mdata) > 0:
            pd = dict([(sid, dd['CA']) for sid, dd in mdata.iteritems()])
            self.add_section(
                name = 'CpA Retention in Each Read',
                anchor = 'biscuit-retention',
                description = "<p>This plot shows the distribution of the number of retained CpA cytosine in each read, up to 10.</p>",
                plot = linegraph.plot(pd, {'ylab': 'Number of Reads', 'xlab': 'Number of Retention within Read'})
            )

            pd = dict([(sid, dd['CC']) for sid, dd in mdata.iteritems()])
            self.add_section(
                name = 'CpC Retention in Each Read',
                anchor = 'biscuit-retention',
                description = "<p>This plot shows the distribution of the number of retained CpC cytosine in each read, up to 10.</p>",
                plot = linegraph.plot(pd, {'ylab': 'Number of Reads', 'xlab': 'Number of Retention within Read'})
            )

            pd = dict([(sid, dd['CG']) for sid, dd in mdata.iteritems()])
            self.add_section(
                name = 'CpG Retention in Each Read',
                anchor = 'biscuit-retention',
                description = "<p>This plot shows the distribution of the number of retained CpG cytosine in each read, up to 10.</p>",
                plot = linegraph.plot(pd, {'ylab': 'Number of Reads', 'xlab': 'Number of Retention within Read'})
            )

            pd = dict([(sid, dd['CT']) for sid, dd in mdata.iteritems()])
            self.add_section(
                name = 'CpT Retention in Each Read',
                anchor = 'biscuit-retention',
                description = "<p>This plot shows the distribution of the number of retained CpT cytosine in each read, up to 10.</p>",
                plot = linegraph.plot(pd, {'ylab': 'Number of Reads', 'xlab': 'Number of Retention within Read'})
            )

        mdata = dict([(k.replace('_CpHRetentionByReadPos.txt',''),v['1']) for k, v in self.mdata['retention'].iteritems() if k.endswith('_CpHRetentionByReadPos.txt')])
        if len(mdata) > 0 and all([len(v)>0 for v in mdata.itervalues()]):
            self.add_section(
                name = 'CpH Retention by Position in Read 1',
                anchor = 'biscuit-retention',
                description = "<p>This plot shows the distribution of CpH retention rate in read 1.</p>",
                plot = linegraph.plot(mdata, {'ylab': 'CpH Retention Rate (%)', 'xlab': 'Position in Read'})
            )

        mdata = dict([(k.replace('_CpHRetentionByReadPos.txt',''),v['2']) for k, v in self.mdata['retention'].iteritems() if k.endswith('_CpHRetentionByReadPos.txt')])
        if len(mdata) > 0 and all([len(v)>0 for v in mdata.itervalues()]):
            self.add_section(
                name = 'CpH Retention by Position in Read 2',
                anchor = 'biscuit-retention',
                description = "<p>This plot shows the distribution of CpH retention rate in read 2.</p>",
                plot = linegraph.plot(mdata, {'ylab': 'CpH Retention Rate (%)', 'xlab': 'Position in Read'})
            )
