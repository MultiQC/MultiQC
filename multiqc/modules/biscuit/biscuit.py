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

def _uniquify_ID(data, header, namespace):
    header2 = dict([(namespace+k, v) for k, v in header.items()])
    data2 = {}
    for sid, datum in data.items():
        data2[sid] = dict([(namespace+k,v) for k, v in datum.items()])

    return (data2, header2)

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
                if k != 'markdup':
                    getattr(self, 'chart_biscuit_%s' % k)()

    def biscuit_stats_table(self):

        pd = {}
        for sid, dd in self.mdata['mapq'].items():
            if not sid.endswith('_mapq_table.txt'):
                continue
            sid = sid.replace('_mapq_table.txt','')
            allreads = sum([int(_) for _ in dd.values()])
            pd[sid] = {'%aligned': float(allreads-int(dd['unmapped']))/allreads*100}
        
        hasPE = False
        hasSE = False
        for sid, dd in self.mdata['markdup'].items():
            sid = sid.replace('_markdup_report.txt','')
            sid = sid.replace('_markdup.bam','')
            sid = sid.replace('.bam','')
            sid = sid.replace('_markdup','')
            if sid not in pd:
                pd[sid] = {}
            if 'dupRatePE' in dd and dd['dupRatePE'] is not None:
                pd[sid]['%dupratePE'] = dd['dupRatePE']
                hasPE = True
            if 'dupRateSE' in dd and dd['dupRateSE'] is not None:
                pd[sid]['%duprateSE'] = dd['dupRateSE']
                hasSE = True

        pheader = OrderedDict()
        pheader['%aligned'] = {'title':'% Aligned', 'max':100, 'min':0, 'suffix':'%','scale':'Greens'}
        if hasPE:
            pheader['%dupratePE'] = {'title':'% DupsPE', 'max':100, 'min':0, 'suffix':'%','scale':'Reds'}
        if hasSE:
            pheader['%duprateSE'] = {'title':'% DupsSE', 'max':100, 'min':0, 'suffix':'%','scale':'Reds'}
        self.general_stats_addcols(pd, pheader)

        # if len(mdata)>0 and all([dd['Duplication Rate'] is not None for dd in mdata.values()]):
        # self.general_stats_addcols(mdata, {'% Dup': {'suffix':'%', 'min':0, 'max':100}})


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

        if fn.endswith('_isize_score_table.txt'):
            data['I'] = {}      # insert size
            data['S'] = {}      # SW-score

            for l in f.splitlines():
                fields = l.split()
                if fields[0] == 'I':
                    data[fields[0]][int(fields[1])] = float(fields[2]) * 100
                elif fields[0] == 'S':
                    data[fields[0]][fields[1]] = float(fields[2])
                
        return data
    
    def chart_biscuit_mapq(self):

        # fraction of optimally mapped reads
        pd = {}
        for sid, dd in self.mdata['mapq'].items():
            if not sid.endswith('_mapq_table.txt'):
                continue
            sid = sid.replace('_mapq_table.txt','');
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
            plot = bargraph.plot(pd, OrderedDict([
                ('OAligned', {'name':'Optimally Aligned Reads', 'color': '#8bbc21'}),
                ('SAligned', {'name':'Suboptimally Aligned Reads', 'color': '#f7a35c'}),
                ('UAligned', {'name':'Unaligned Reads', 'color': '#000000'})
            ]), {'id':'biscuit_mapping',
                 'title':'BISCUIT: Mapping Summary',
                 'ylab':'Number of Reads',
                 'cpswitch_c_active': True,
                 'cpswitch_counts_label': '# Reads'
            })
        )

        # mapping strand distribution
        pd = dict([(k.replace('_strand_table.txt', ''),v) for k,v in self.mdata['mapq'].items() if k.endswith('_strand_table.txt')])
        self.add_section(
            name='Mapping Strand Distribution',
            anchor='biscuit-strands',
            description = "<p>This plot shows the distribution of strand of mapping and strand of bisulfite conversion.</p>",
            plot = bargraph.plot(pd, OrderedDict([
                ('++', {'name':'Waston-Aligned, Waston-Bisulfite Conversion', 'color': '#F53855'}),
                ('+-', {'name':'Waston-Aligned, Crick-Bisulfite Conversion', 'color': '#E37B40'}),
                ('-+', {'name':'Crick-Aligned, Waston-Bisulfite Conversion', 'color': '#46B29D'}),
                ('--', {'name':'Crick-Aligned, Crick-Bisulfite Conversion', 'color': '#324D5C'})
            ]), {'id':'biscuit_strands',
                 'title':'BISCUIT: Mapping Strand Distribution',
                 'ylab':'Number of Reads',
                 'cpswitch_c_active': True,
                 'cpswitch_counts_label': '# Reads'
            })
        )

        # Mapping Quality, Insert Size
        total = {}
        for sid, dd in self.mdata['mapq'].items():
            if not sid.endswith('_mapq_table.txt'):
                continue
            sid = sid.replace('_mapq_table.txt','')
            total[sid] = sum([int(cnt) for _, cnt in dd.items() if _ != "unmapped"])

        pd_mapping = {}
        for sid, dd in self.mdata['mapq'].items():
            if not sid.endswith('_mapq_table.txt'):
                continue
            sid = sid.replace('_mapq_table.txt','')
            mapqcnts = []
            for mapq in range(61):
                if str(mapq) in dd:
                    mapqcnts.append(float(dd[str(mapq)])/total[sid]*100)
                else:
                    mapqcnts.append(0)
            pd_mapping[sid] = dict(zip(range(61), mapqcnts))

        pd_isize = dict([(k.replace('_isize_score_table.txt', ''),v['I']) for k,v in self.mdata['mapq'].items() if k.endswith('_isize_score_table.txt')])

        self.add_section(
            name = 'Mapping Quality and Insert Size',
            anchor = 'biscuit-mapq-isize',
            description = "<p>This plot shows the distribution of Primary Mapping Quality.</p>",
            plot = linegraph.plot([pd_mapping, pd_isize],
                {'id':'biscuit_mapping', 'title': 'Mapping Information', 'data_labels': [
                    {'name':'Mapping Quality', 'ylab': '% Primary Mapped Reads','xlab': 'Mapping Quality'},
                    {'name':'Insert Size', 'ylab': '% Mapped Reads', 'xlab': 'Insert Size'}]}))

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
        mdata = [
            dict([(k.replace('_covdist_table.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_covdist_table.txt')]),
            dict([(k.replace('_covdist_q40_table.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_covdist_q40_table.txt')]),
            dict([(k.replace('_covdist_q40_botgc_table.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_covdist_q40_botgc_table.txt')]),
            dict([(k.replace('_covdist_q40_topgc_table.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_covdist_q40_topgc_table.txt')]),
            dict([(k.replace('_covdist_cpg_table.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_covdist_cpg_table.txt')]),
            dict([(k.replace('_covdist_cpg_q40_table.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_covdist_cpg_q40_table.txt')]),
            dict([(k.replace('_covdist_cpg_q40_botgc_table.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_covdist_cpg_q40_botgc_table.txt')]),
            dict([(k.replace('_covdist_cpg_q40_topgc_table.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_covdist_cpg_q40_topgc_table.txt')]),
        ]
        if len(mdata)>0:
            self.add_section(
                name = 'Cumulative Base Coverage',
                anchor = 'biscuit-coverage-base',
                description = "<p>This plot shows the cummulative base coverage. High and low GC content region are the top and bottom 10% 100bp window in GC content.</p>",
                plot = linegraph.plot(mdata, {'id':'biscuit_coverage_base',
                    'title': 'Cumulative Base Coverage',
                    'data_labels': [
                        {'name': 'All', 'ylab':'Million Bases'}, 
                        {'name': 'Q40', 'ylab':'Million Bases'}, 
                        {'name': 'Q40 low GC', 'ylab':'Million Bases'}, 
                        {'name': 'Q40 high GC', 'ylab':'Million Bases'}, 
                        {'name': 'CpG (all)', 'ylab':'Million CpGs'}, 
                        {'name': 'CpG Q40', 'ylab':'Million CpGs'}, 
                        {'name': 'CpG Q40 low GC', 'ylab':'Million CpGs'}, 
                        {'name': 'CpG Q40 high GC', 'ylab':'Million CpGs'}, 
                    ]})
            )
            
            for sid, dd in mdata[0].items():
                if sid not in basecov:
                    basecov[sid] = {}
                basecov[sid]['all'] = dd[1]/dd[0]*100

            for sid, dd in mdata[1].items():
                if sid not in basecov:
                    basecov[sid] = {}
                basecov[sid]['uniq'] = dd[1]/dd[0]*100

            for sid, dd in mdata[4].items():
                if sid not in basecov:
                    basecov[sid] = {}
                basecov[sid]['cpg_all'] = dd[1]/dd[0]*100            

            for sid, dd in mdata[5].items():
                if sid not in basecov:
                    basecov[sid] = {}
                basecov[sid]['cpg_uniq'] = dd[1]/dd[0]*100

        # base coverage >=1x table
        if len(basecov)>0:
            self.add_section(
                name = 'Coverage by At Least One Read',
                anchor = 'biscuit-coverage-base-table',
                description = '<p>The fraction covered by at least one read.</p>',
                plot = table.plot(basecov, {
                    'all':{'title':'Genome (All)','max':100,'min':0,'suffix':'%'},
                    'uniq':{'title':'Genome (Unique)','max':100,'min':0,'suffix':'%'},
                    'cpg_all':{'title':'CpG (All)','max':100,'min':0,'suffix':'%'},
                    'cpg_uniq':{'title':'CpG (Unique)','max':100,'min':0,'suffix':'%'},
                })
            )

        # cpg distribution
        mdata = dict([(k.replace('_cpg_dist_table.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_cpg_dist_table.txt')])
        if len(mdata)>0:
            hdr = OrderedDict()
            pd = OrderedDict()
            sid = 'Genome'
            dd = list(mdata.values())[0]
            pd[sid] = OrderedDict()
            for ctg in ['ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
                ctg1 = ctg.replace('CpGs','')
                hdr[ctg1] = {'max':100,'min':0,'suffix':'%'}
                pd[sid][ctg1] = float(dd[ctg]['uc']) / dd['TotalCpGs']['uc'] * 100

            hdr['Exonic']['description'] = 'Exonic CpGs'
            hdr['Repeat']['description'] = 'Repeat-Masked CpGs'
            hdr['Genic']['description']  = 'Genic CpGs'
            hdr['CGI']['description']    = 'CpG Island CpGs'
            for sid, dd in mdata.items():
                pd[sid] = OrderedDict()
                for ctg in ['ExonicCpGs', 'RepeatCpGs', 'GenicCpGs', 'CGICpGs']:
                    ctg1 = ctg.replace('CpGs','')
                    pd[sid][ctg1] = float(dd[ctg]['uc']) / dd['TotalCpGs']['uc'] * 100

            self.add_section(
                name = 'CpG Coverage vs Genomic Features',
                anchor = 'biscuit-coverage-cpg-dist',
                description = "<p>The top row shows how CpGs breaks down to different categories. Each other row shows the how CpGs uniquely covered by the given data breaks down to these categories. It is the fraction of CpGs in the given category out of all CpGs covered by the data.</p>",
                plot = table.plot(pd, hdr)
            )

            pd = dict([(sid, dd['cgi_coverage']) for sid, dd in mdata.items()])
            self.add_section(
                name = 'CpG Island Coverage',
                anchor = 'biscuit-coverage-cgi',
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
                ]), {'id':'cgi-cov-table'})
            )
        
        ####################
        # sequencing depth
        ####################
        pd = OrderedDict()
        mdata = dict([(k.replace('_all_cv_table.txt',''),v) 
                for k, v in self.mdata['coverage'].items() if k.endswith('_all_cv_table.txt')])
        for sid, dd in mdata.items():
            pd[sid] = OrderedDict()
            for ctg in ['all','all_topgc','all_botgc']:
                if ctg in dd:
                    pd[sid]['cv_'+ctg] = dd[ctg]['cv']
                    pd[sid]['mu_'+ctg] = dd[ctg]['mu']
                    #  pd[sid]['sigma_'+ctg] = dd[ctg]['sigma']

        mdata = dict([(k.replace('_cpg_cv_table.txt',''),v) 
                for k, v in self.mdata['coverage'].items() if k.endswith('_cpg_cv_table.txt')])
        for sid, dd in mdata.items():
            if sid not in pd:
                pd[sid] = OrderedDict()
            for ctg in ['cpg','cpg_topgc','cpg_botgc']:
                if ctg in dd:
                    pd[sid]['cv_'+ctg] = dd[ctg]['cv']
                    pd[sid]['mu_'+ctg] = dd[ctg]['mu']
                    #  pd[sid]['sigma_'+ctg] = dd[ctg]['sigma']

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
        #  pheader['sigma_all'] = {'title':'Genome','description':'Whole Genome Average'}
        #  pheader['sigma_all_topgc'] = {'title':'Highest GC','description':'Top Decile in GC Content'}
        #  pheader['sigma_all_botgc'] = {'title':'Lowest GC','description':'Bottom Decile in GC Content'}

        self.add_section(
            name = 'Sequencing Depth',
            anchor = 'biscuit-seq-depth',
            description = "<p>This plot shows sequence depth mean and uniformity measured in Coefficient of Variation (mu/sigma), mapQ>40 only. CG.* shows the depth on CpG only. GC contents were measured on 100bp non-overlapping windows.</p>",
            plot = table.plot(pd, pheader))

        




        ############################################################
        ### duplicate read coverage
        ### this is commented out to reduce the number of sections 
        #############################################################
        #  mdata = dict([(k.replace('_dup_report.txt',''),v) for k, v in self.mdata['coverage'].items() if k.endswith('_dup_report.txt')])
        #  self.add_section(
        #      name = 'Duplicate Read Coverage',
        #      anchor = 'biscuit-coverage-dup',
        #      description = "<p>This table shows base coverage by duplicate reads as a percentage of coverage by all reads.</p>",
        #      plot = table.plot(*_uniquify_ID(mdata, OrderedDict([
        #          ('all', {'id':'cov-dup-all', 'title':'Genome', 'suffix':'%', 'description':'Fraction of coverage by duplicate reads'}),
        #          ('topGC', {'id':'cov-dup-topgc', 'title':'HighGC', 'suffix':'%', 'description':'Fraction of coverage in top decile in GC content by duplicate reads'}),
        #          ('botGC', {'id':'cov-dup-botgc', 'title':'LowGC', 'suffix':'%', 'description':'Fraction of coverage in bottom decile in GC content by duplicate reads'}),
        #          ('all-q40', {'id':'cov-all-q40', 'title':'Genome-uniq', 'suffix':'%', 'description':'Fraction of coverage by duplicate reads with mapq>=40'}),
        #          ('topGC-q40', {'id':'cov-dup-topgc-q40', 'title':'HighGC-uniq', 'suffix':'%',
        #           'description':'Fraction of coverage in top decile in GC content by duplicate reads with mapq>=40'}),
        #          ('botGC-q40', {'id':'cov-dup-botgc-q40', 'title':'LowGC-uniq', 'suffix':'%',
        #           'description':'Fraction of coverage in bottom decile in GC content by duplicate reads with mapq>=40'}),
        #      ]), 'biscuit-coverage-dup')))

    def parse_logs_retention(self, f, fn):

        data = {}
        if fn.endswith('_totalReadConversionRate.txt'):
            m = re.match(r'([\d.]+)\t([\d.]+)\t([\d.]+)\t([\d.]+)', f.splitlines()[1])
            data['ca'] = float(m.group(1))*100
            data['cc'] = float(m.group(2))*100
            data['cg'] = float(m.group(3))*100
            data['ct'] = float(m.group(4))*100

        elif fn.endswith('_totalBaseConversionRate.txt'):
            m = re.match(r'([\d.]+)\t([\d.]+)\t([\d.]+)\t([\d.]+)', f.splitlines()[1])
            data['bca'] = float(m.group(1))*100
            data['bcc'] = float(m.group(2))*100
            data['bcg'] = float(m.group(3))*100
            data['bct'] = float(m.group(4))*100

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

        elif fn.endswith('_CpHRetentionByReadPos.txt') or fn.endswith('_CpGRetentionByReadPos.txt'):
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

        elif fn.endswith('_beta_table.txt'):
            sumcounts = 0
            for l in f.splitlines():
                fields = l.split('\t')
                sumcounts += int(fields[1])
            
            data = {}
            for l in f.splitlines():
                fields = l.split('\t')
                data[int(fields[0])] = int(fields[1]) / float(sumcounts)
            
        return data

    def chart_biscuit_retention(self):

        mdata1 = dict([(k.replace('_totalReadConversionRate.txt',''),v) for k, v in self.mdata['retention'].items() if k.endswith('_totalReadConversionRate.txt')])
        mdata2 = dict([(k.replace('_totalBaseConversionRate.txt',''),v) for k, v in self.mdata['retention'].items() if k.endswith('_totalBaseConversionRate.txt')])
        pdata = dict([(sid, dict(list(dd.items()) + list(mdata2[sid].items()))) for sid, dd in mdata1.items()])

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
            description = "<p>This plot shows cytosine retention rate. 'r.' stands for read-averaging and 'b.' stands for 'base-averaging.</p>",
            plot = table.plot(pdata, pheader))

        ## number of cytosine retention distribution
        mdata = dict([(k.replace('_freqOfTotalRetentionPerRead.txt',''),v) for k, v in self.mdata['retention'].items() if k.endswith('_freqOfTotalRetentionPerRead.txt')])
        mdata_meth = dict([(k.replace('_beta_table.txt',''),v) 
                for k,v in self.mdata['retention'].items() if k.endswith('_beta_table.txt')])
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
            description = "<p>This plot shows the distribution of the number of retained cytosine in each read, up to 10.</p>",
            plot = linegraph.plot(pd, {
                'id': 'biscuit_retention_read_cpa', 
                'xlab': 'Number of Retention within Read',
                'title': 'Retention Distribution',
                'data_labels': [
                    {'name': 'CpG retention', 'ylab': 'Fraction of cytosine in CpG context', 'xlab': 'Retention Level (%)'},
                    {'name': 'Within-read CpA', 'ylab': 'Number of Reads'},
                    {'name': 'Within-read CpC', 'ylab': 'Number of Reads'},
                    {'name': 'Within-read CpG', 'ylab': 'Number of Reads'},
                    {'name': 'Within-read CpT', 'ylab': 'Number of Reads'},
                ]})
            )

        ## retention vs read position
        mdata = [
            dict([(k.replace('_CpHRetentionByReadPos.txt',''),v['1']) for k, v in self.mdata['retention'].items() if k.endswith('_CpHRetentionByReadPos.txt')]),
            dict([(k.replace('_CpHRetentionByReadPos.txt',''),v['2']) for k, v in self.mdata['retention'].items() if k.endswith('_CpHRetentionByReadPos.txt')]),
            dict([(k.replace('_CpGRetentionByReadPos.txt',''),v['1']) for k, v in self.mdata['retention'].items() if k.endswith('_CpGRetentionByReadPos.txt')]),
            dict([(k.replace('_CpGRetentionByReadPos.txt',''),v['2']) for k, v in self.mdata['retention'].items() if k.endswith('_CpGRetentionByReadPos.txt')]),
        ]
        self.add_section(
            name = 'Retention vs. Base Position in Read',
            anchor = 'biscuit-retention-cytosine',
            description = "<p>This plot (aka. mbias plot) shows the distribution of cytosine retention rate in read.</p>",
            plot = linegraph.plot(mdata, {
                'id': 'biscuit_retention_cytosine',
                'xlab': 'Position in Read', 'ymin':0, 'ymax':100, 'yMinRange':0, 'yFloor':0,
                'title': 'Retention vs. Base Position in Read',
                'data_labels': [
                    {'name': 'CpH Read 1', 'ylab': 'CpH Retention Rate (%)', 'ymin':0, 'ymax':100},
                    {'name': 'CpH Read 2', 'ylab': 'CpH Retention Rate (%)', 'ymin':0, 'ymax':100},
                    {'name': 'CpG Read 1', 'ylab': 'CpG Retention Rate (%)', 'ymin':0, 'ymax':100},
                    {'name': 'CpG Read 2', 'ylab': 'CpG Retention Rate (%)', 'ymin':0, 'ymax':100},
                ]})
            )

