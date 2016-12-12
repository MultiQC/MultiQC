#!/usr/bin/env python

""" MultiQC module to parse output from Slamdunk """

from __future__ import print_function
import logging
import re
from distutils.version import StrictVersion
from collections import OrderedDict

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Slamdunk module class, parses slamdunk logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Slamdunk', anchor='slamdunk',
        href='http://t-neumann.github.io/slamdunk/',
        info="is a tool to analyze SLAMSeq data.")

        self.slamdunk_data = dict()
        
        self.PCA_data = dict()
        
        self.utrates_data = dict()
        
        self.rates_data_plus = dict()
        self.rates_data_minus = dict()
        
        self.nontc_per_readpos_plus = dict()
        self.nontc_per_readpos_minus = dict()
        
        self.tc_per_readpos_plus = dict()
        self.tc_per_readpos_minus = dict()
        
        self.nontc_per_utrpos_plus = dict()
        self.nontc_per_utrpos_minus = dict()
        
        self.tc_per_utrpos_plus = dict()
        self.tc_per_utrpos_minus = dict()

        # List of plots to produce
        produce = []
            
        for f in self.find_log_files(config.sp['slamdunk']['summary'], filehandles = True):
            self.parseSummary(f)
            
        if len(self.slamdunk_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
        else :
            produce.append("summary")
            
        log.debug("Parsing PCA reports.")
            
        for f in self.find_log_files(config.sp['slamdunk']['PCA'], filehandles = True):
            self.parsePCA(f)
            
        if len(self.PCA_data) == 0:
            log.debug("Could not find PCA reports in {}".format(config.analysis_dir))
        else :
            produce.append("PCA")
            
        log.debug("Parsing UTR rates reports. This can take a while...")
                        
        for f in self.find_log_files(config.sp['slamdunk']['utrrates'], filehandles = True):
            self.parseUtrRates(f)
            
        if len(self.utrates_data) == 0:
            log.debug("Could not find UTR rates reports in {}".format(config.analysis_dir))
        else :
            self.write_data_file(self.utrates_data, 'multiqc_slamdunk_utrrates')
            produce.append("utrates")            

        log.debug("Parsing read rates reports.")
            
        for f in self.find_log_files(config.sp['slamdunk']['rates'], filehandles = True):
            self.parseSlamdunkRates(f)
            
        if len(self.rates_data_plus) == 0:
            log.debug("Could not find read rates reports in {}".format(config.analysis_dir))
        else :
            self.write_data_file(self.rates_data_plus, 'multiqc_slamdunk_readrates_plus')
            self.write_data_file(self.rates_data_minus, 'multiqc_slamdunk_readrates_minus')
            produce.append("readrates")
            
        log.debug("Parsing rates per read position reports.")
            
        for f in self.find_log_files(config.sp['slamdunk']['tcperreadpos'], filehandles = True):
            self.parseSlamdunkTCPerReadpos(f)
            
        if len(self.tc_per_readpos_plus) == 0:
            log.debug("Could not find conversion per read position reports in {}".format(config.analysis_dir))
        else :
            self.write_data_file(self.tc_per_readpos_plus, 'multiqc_slamdunk_tcperreadpos_plus')
            self.write_data_file(self.nontc_per_readpos_plus, 'multiqc_slamdunk_nontcperreadpos_plus')
            self.write_data_file(self.tc_per_readpos_minus, 'multiqc_slamdunk_tcperreadpos_minus')
            self.write_data_file(self.nontc_per_readpos_minus, 'multiqc_slamdunk_nontcperreadpos_minus')
            produce.append("readpos")
            
        log.debug("Parsing rates per UTR position reports.")

        for f in self.find_log_files(config.sp['slamdunk']['tcperutrpos'], filehandles = True):
            self.parseSlamdunkTCPerUtrpos(f)
            
        if len(self.nontc_per_utrpos_plus) == 0:
            log.debug("Could not find conversion per UTR position reports in {}".format(config.analysis_dir))
        else :
            self.write_data_file(self.tc_per_utrpos_plus, 'multiqc_slamdunk_tcperutrpos_plus')
            self.write_data_file(self.nontc_per_utrpos_plus, 'multiqc_slamdunk_nontcperutrpos_plus')
            self.write_data_file(self.tc_per_utrpos_minus, 'multiqc_slamdunk_tcperutrpos_minus')
            self.write_data_file(self.nontc_per_utrpos_minus, 'multiqc_slamdunk_nontcperutrpos_minus')
            produce.append("utrpos")
        
        if not produce:
            log.debug("No slamdunk reports found.")
            raise UserWarning
        
        # Start the sections
        self.sections = list()
        
        # Basic Stats Table
        if "summary" in produce:
            self.slamdunkGeneralStatsTable()
            self.slamdunkFilterStatsTable()
        
        # PCA plot
        if "PCA" in produce:
            self.slamdunkPCAPlot()
        
        # Utr rates plot
        if "utrates" in produce:
            self.slamdunkUtrRatesPlot()

        # Rates plot
        if "readrates" in produce:
            self.slamdunkOverallRatesPlot()
        
        # TC per read position plot
        if "readpos" in produce:
            self.slamdunkTcPerReadPosPlot()
        
        # TC per UTR position plot
        if "utrpos" in produce:
            self.slamdunkTcPerUTRPosPlot()
            
        log.info("Found " + str(len(self.slamdunk_data)) + " reports")
        
    def parsePCA(self, f):
        
        # Skip header
        next(f['f'])
        
        for line in f['f']:
            fields = line.rstrip().split('\t')
            
            sample = fields[0]
            PC1 = fields[1]
            PC2 = fields[2]
            
            self.PCA_data[sample] = dict()

            self.PCA_data[sample] = [{'x': float(PC1), 'y': float(PC2)}]
        
        
    def parseUtrRates(self, f) :
                
        # Skip comment line #
        next(f['f'])
        
        # Read median header
        line = next(f['f'])
        
        if "Conversions=" in line :
            
            sample = f['s_name']
            self.utrates_data[sample] = OrderedDict()
            
            conversions = re.sub(".*Conversions=","",line.rstrip()).split(",")
            
            for conversion in conversions:
                type, value = conversion.split(":")
                self.utrates_data[sample][type] = float(value)
        
        else :
            
            log.warning("Malformed UTR rates header. Conversion rates per UTR plot will be affected.")

    def parseSlamdunkRates(self, f):
        
        sample = f['s_name']
        
        # Skip comment line #
        next(f['f'])
        
        bases = next(f['f']).rstrip().split('\t')
        
        baseDict = {}
        order = {}
        
        for i in range(1, len(bases)) :
            order[i] = bases[i]
                    
        for line in f['f']:
            values = line.rstrip().split('\t')
            base = values[0]
            baseDict[base]= {}
            
            for i in range(1, len(values)) :
            
                baseDict[base][order[i]] = int(values[i])
                
        divisor = {}      
        
        for fromBase in baseDict:
            for toBase in baseDict[fromBase]:
                if(toBase.islower()) :
                    if not divisor.has_key(fromBase.lower()) :
                        divisor[fromBase.lower()] = 0
                    divisor[fromBase.lower()] += baseDict[fromBase][toBase]
                else:
                    if not divisor.has_key(fromBase) :
                        divisor[fromBase] = 0
                    divisor[fromBase] += baseDict[fromBase][toBase]
        
        for fromBase in baseDict:
            for toBase in baseDict[fromBase]:
                if(toBase.islower()) :
                    baseDict[fromBase][toBase] = baseDict[fromBase][toBase] / float(divisor[fromBase.lower()]) * 100
                else:
                    baseDict[fromBase][toBase] = baseDict[fromBase][toBase] / float(divisor[fromBase]) * 100
        
        self.rates_data_plus[sample] = {}
        self.rates_data_minus[sample] = {}
        
        for fromBase in baseDict:
            for toBase in baseDict[fromBase]:
                if fromBase != "N" and toBase.upper() != "N" and fromBase != toBase.upper():
                    if(toBase.islower()) :
                        self.rates_data_minus[sample][fromBase + ">" + toBase.upper()] = baseDict[fromBase][toBase]
                    else :
                        self.rates_data_plus[sample][fromBase + ">" + toBase] = baseDict[fromBase][toBase]
                        
    def parseSlamdunkTCPerReadpos(self, f):
        
        sample = f['s_name']
        
        # Skip comment line #
        next(f['f'])
        
        self.nontc_per_readpos_plus[sample] = {}
        self.nontc_per_readpos_minus[sample] = {}
        
        self.tc_per_readpos_plus[sample] = {}
        self.tc_per_readpos_minus[sample] = {}
        
        pos = 1
        
        for line in f['f']:
            values = line.rstrip().split('\t')
            if int(values[4]) > 0 :
                self.nontc_per_readpos_plus[sample][pos] = float(int(values[0])) / int(values[4]) * 100
                self.tc_per_readpos_plus[sample][pos] = float(int(values[2])) / int(values[4]) * 100
            else :
                self.nontc_per_readpos_plus[sample][pos] = 0
                self.tc_per_readpos_plus[sample][pos] = 0
                
            if int(values[5]) > 0:
                self.nontc_per_readpos_minus[sample][pos] = float(int(values[1])) / int(values[5]) * 100
                self.tc_per_readpos_minus[sample][pos] = float(int(values[3])) / int(values[5]) * 100
            else:
                self.nontc_per_readpos_minus[sample][pos] = 0
                self.tc_per_readpos_minus[sample][pos] = 0
                
            pos += 1
            
    def parseSlamdunkTCPerUtrpos(self, f):
        
        sample = f['s_name']
        
        # Skip comment line #
        next(f['f'])
        
        self.nontc_per_utrpos_plus[sample] = {}
        self.nontc_per_utrpos_minus[sample] = {}
        
        self.tc_per_utrpos_plus[sample] = {}
        self.tc_per_utrpos_minus[sample] = {}
        
        pos = 1
        
        for line in f['f']:
            values = line.rstrip().split('\t')
            if int(values[4]) > 0 :
                self.nontc_per_utrpos_plus[sample][pos] = float(int(values[0])) / int(values[4]) * 100
                self.tc_per_utrpos_plus[sample][pos] = float(int(values[2])) / int(values[4]) * 100
            else :
                self.nontc_per_utrpos_plus[sample][pos] = 0
                self.tc_per_utrpos_plus[sample][pos] = 0
                
            if int(values[5]) > 0:
                self.nontc_per_utrpos_minus[sample][pos] = float(int(values[1])) / int(values[5]) * 100
                self.tc_per_utrpos_minus[sample][pos] = float(int(values[3])) / int(values[5]) * 100
            else:
                self.nontc_per_utrpos_minus[sample][pos] = 0
                self.tc_per_utrpos_minus[sample][pos] = 0
                
            pos += 1
    
    def parseSummary(self, f):
        
        # Skip comment line #
        next(f['f'])
        
        # Skip header line "FileName..."
        columnCount = next(f['f']).count("\t") + 1
            
        for line in f['f']:

            fields = line.rstrip().split("\t")
            self.slamdunk_data[self.clean_s_name(fields[0],"")] = dict()
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['sequenced'] = int(fields[4])
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['mapped'] = int(fields[5])
            #self.slamdunk_data[self.clean_s_name(fields[0],"")]['deduplicated'] = int(fields[6])
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['mqfiltered'] = int(fields[7])
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['idfiltered'] = int(fields[8])
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['nmfiltered'] = int(fields[9])
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['multimapper'] = int(fields[10])
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['retained'] = int(fields[11])
            
            # Additional Count Column found in Table
            if columnCount == 14:
                self.slamdunk_data[self.clean_s_name(fields[0],"")]['counted'] = int(fields[12])
                
        self.add_data_source(f)
        
    def slamdunkGeneralStatsTable(self):
        """ Take the parsed summary stats from Slamdunk and add it to the
        basic stats table at the top of the report """


        headers = OrderedDict()
        
        headers['counted'] = {
            'title': 'Counted',
            'description': '# reads counted within 3\'UTRs',
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:.f}',
            'scale': 'YlGn',
            'modify': lambda x: x / 1000000,
        }
        
        headers['retained'] = {
            'title': 'Retained',
            'description': '# retained reads after filtering',
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:.f}',
            'scale': 'YlGn',
            'modify': lambda x: x / 1000000,
        }
        
#         headers['multimapper'] = {
#             'title': 'Multimap-Filtered',
#             'description': '# multimap-filtered reads',
#             'shared_key': 'read_count',
#             'min': 0,
#             'format': '{:.f}',
#             'scale': 'OrRd',
#         }
#         
#         headers['nmfiltered'] = {
#             'title': 'NM-Filtered',
#             'description': '# NM-filtered reads',
#             'shared_key': 'read_count',
#             'min': 0,
#             'format': '{:.f}',
#             'scale': 'OrRd',
#         }
#         
#         headers['idfiltered'] = {
#             'title': 'Identity-Filtered',
#             'description': '# identity-filtered reads',
#             'shared_key': 'read_count',
#             'min': 0,
#             'format': '{:.f}',
#             'scale': 'OrRd',
#         }
#         
#         headers['mqfiltered'] = {
#             'title': 'MQ-Filtered',
#             'description': '# MQ-filtered reads',
#             'shared_key': 'read_count',
#             'min': 0,
#             'format': '{:.f}',
#             'scale': 'OrRd',
#         }
        
        headers['mapped'] = {
            'title': 'Mapped',
            'description': '# mapped reads',
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:.f}',
            'scale': 'YlGn',
            'modify': lambda x: x / 1000000,
        }
        
        headers['sequenced'] = {
            'title': 'Sequenced',
            'description': '# sequenced reads',
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:.f}',
            'scale': 'YlGn',
            'modify': lambda x: x / 1000000,
        }
        
        self.general_stats_addcols(self.slamdunk_data, headers)
        
    def slamdunkFilterStatsTable(self):
        """ Take the parsed filter stats from Slamdunk and add it to a separate table """

        headers = OrderedDict()
        
        headers['mapped'] = {
            'title': 'Mapped',
            'description': '# mapped reads',
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:.f}',
            'scale': 'YlGn',
            'modify': lambda x: x / 1000000,
        }
        
        headers['multimapper'] = {
            'title': 'Multimap-Filtered',
            'description': '# multimap-filtered reads',
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:.f}',
            'scale': 'OrRd',
            'modify': lambda x: x / 1000000,
        }
        
        headers['nmfiltered'] = {
            'title': 'NM-Filtered',
            'description': '# NM-filtered reads',
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:.f}',
            'scale': 'OrRd',
            'modify': lambda x: x / 1000000,
        }
        
        headers['idfiltered'] = {
            'title': 'Identity-Filtered',
            'description': '# identity-filtered reads',
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:.f}',
            'scale': 'OrRd',
            'modify': lambda x: x / 1000000,
        }
        
        headers['mqfiltered'] = {
            'title': 'MQ-Filtered',
            'description': '# MQ-filtered reads',
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:.f}',
            'scale': 'OrRd',
            'modify': lambda x: x / 1000000,
        }
        
        config = {
            'id': 'slamdunk_filtering_table',
            'min': 0,
        }
        
        self.sections.append({
            'name': 'Filter statistics',
            'anchor': 'slamdunk_filtering',
            'content': '<p> This table shows the number of reads filtered with each filter criterion during filtering phase of slamdunk.</p>' +  
                        plots.table.plot(self.slamdunk_data, headers, config)
        })
                
    def slamdunkOverallRatesPlot (self):
        """ Generate the overall rates plot """

        pconfig = {
            'id': 'overallratesplot',
            'title': 'Overall conversion rates in reads',
            'cpswitch' : False,
            'cpswitch_c_active': False,
            'stacking' : 'normal',
            'tt_decimals': 2,
            'tt_suffix': '%',
            'tt_percentages': False,
            'data_labels': [
                "Plus Strand +",
                "Minus Strand -", 
            ]
        }
        
        cats = [OrderedDict(), OrderedDict()]
        
        cats[0]['T>C'] = {
            #'color': '#D7301F'
            'color': '#fdbf6f'
        }
        cats[0]['A>T'] = {
            #'color': '#C6DBEF'
            'color': '#2171B5'
        }
        cats[0]['A>G'] = {
            #'color': '#6BAED6'
            'color': '#6BAED6'
        }
        cats[0]['A>C'] = {
            #'color': '#2171B5'
            'color': '#C6DBEF'
        }
        cats[0]['T>A'] = {
            #'color': '#C7E9C0'
            'color': '#74C476'
        }
        cats[0]['T>G'] = {
            #'color': '#74C476'
            'color': '#C7E9C0'
        }
        cats[0]['G>A'] = {
            'color': '#D9D9D9'
        }
        cats[0]['G>T'] = {
            'color': '#969696'
        }
        cats[0]['G>C'] = {
            'color': '#525252'
        }
        cats[0]['C>A'] = {
            'color': '#DADAEB'
        }
        cats[0]['C>T'] = {
            'color': '#9E9AC8'
        }
        cats[0]['C>G'] = {
            'color': '#6A51A3'
        }
        
        cats[1]['A>G'] = {
            #'color': '#D7301F'
            'color': '#fdbf6f'
        }
        cats[1]['A>T'] = {
            #'color': '#C6DBEF'
            'color': '#2171B5'
        }
        cats[1]['A>C'] = {
            #'color': '#2171B5'
            'color': '#6BAED6'
        }
        cats[1]['T>A'] = {
            #'color': '#C7E9C0'
            'color': '#C6DBEF'
        }
        cats[1]['T>G'] = {
            #'color': '#74C476'
            'color': '#74C476'
        }
        cats[1]['T>C'] = {
            #'color': '#238B45'
            'color': '#C7E9C0'
        }
        cats[1]['G>A'] = {
            'color': '#D9D9D9'
        }
        cats[1]['G>T'] = {
            'color': '#969696'
        }
        cats[1]['G>C'] = {
            'color': '#525252'
        }
        cats[1]['C>A'] = {
            'color': '#DADAEB'
        }
        cats[1]['C>T'] = {
            'color': '#9E9AC8'
        }
        cats[1]['C>G'] = {
            'color': '#6A51A3'
        }
        
        self.sections.append({
            'name': 'Conversion rates per read',
            'anchor': 'slamdunk_overall_rates',
            'content': '<p>This plot shows the individual conversion rates over all reads. \n\
                        It shows these conversion rates strand-specific: This means for a properly labelled \n\
                        sample you would see a T>C excess on \n\
                        the plus-strand and an A>G excess on the minus strand. <br>\n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#rates" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.bargraph.plot([self.rates_data_plus,self.rates_data_minus], cats, pconfig)
        })
    

    def slamdunkUtrRatesPlot (self):
        """ Generate the UTR rates plot """
        
        cats = OrderedDict()
        cats['T>C'] = {
            #'color': '#D7301F'
            'color': '#fdbf6f'
        }
        cats['A>T'] = {
            #'color': '#C6DBEF'
            'color': '#2171B5'
        }
        cats['A>G'] = {
            #'color': '#6BAED6'
            'color': '#6BAED6'
        }
        cats['A>C'] = {
            #'color': '#2171B5'
            'color': '#C6DBEF'
        }
        cats['T>A'] = {
            #'color': '#C7E9C0'
            'color': '#74C476'
        }
        cats['T>G'] = {
            #'color': '#74C476'
            'color': '#C7E9C0'
        }
        cats['G>A'] = {
            'color': '#D9D9D9'
        }
        cats['G>T'] = {
            'color': '#969696'
        }
        cats['G>C'] = {
            'color': '#525252'
        }
        cats['C>A'] = {
            'color': '#DADAEB'
        }
        cats['C>T'] = {
            'color': '#9E9AC8'
        }
        cats['C>G'] = {
            'color': '#6A51A3'
        }

        pconfig = {
            'id': 'slamdunk_utrratesplot',
            'title': 'Overall conversion rates per UTR',
            'cpswitch' : False,
            'cpswitch_c_active': False,
            'stacking' : 'normal',
            'tt_decimals': 2,
            'tt_suffix': '%',
            'tt_percentages': False,
        }
        
        self.sections.append({
            'name': 'Conversion rates per UTR',
            'anchor': 'slamdunk_utr_rates',
            'content': '<p>This plot shows the individual conversion rates for all UTRs.<br> \n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#utrrates" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.bargraph.plot(self.utrates_data, cats, pconfig)
        })
    
    def slamdunkPCAPlot (self):
        """ Generate the PCA plots """
        
        pconfig = {
            'id': 'slamdunk_pca',
            'title': 'Slamdunk PCA',
            'xlab' : 'PC1',
            'ylab' : 'PC2',
            'tt_label': 'PC1 {point.x:.2f}: PC2 {point.y:.2f}'
        }
        
        self.sections.append({
            'name': 'PCA (T>C based)',
            'anchor': 'slamdunk_PCA',
            'content': '<p> This plot shows the principal components of samples based on the distribution of reads with T>C conversions within UTRs. <br>\n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#summary" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.scatter.plot(self.PCA_data, pconfig) 
        })
        
    def slamdunkTcPerReadPosPlot (self):
        """ Generate the tc per read pos plots """
        
        pconfig_nontc = {
            'id': 'slamdunk_nontcperreadpos_plot',
            'title': 'Non-T>C mutations over reads',
            'ylab': 'Percent mutated %',
            'xlab': 'Position in read',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>Pos {point.x}</b>: {point.y:.2f} %',
            'data_labels': [{'name': 'Forward reads +', 'ylab': 'Percent mutated %'},
                            {'name': 'Reverse reads -', 'ylab': 'Percent mutated %'}]
        }
        
        pconfig_tc = {
            'id': 'slamdunk_tcperreadpos_plot',
            'title': 'T>C conversions over reads',
            'ylab': 'Percent converted %',
            'xlab': 'Position in read',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>Pos {point.x}</b>: {point.y:.2f} %',
            'data_labels': [{'name': 'Forward reads +', 'ylab': 'Percent converted %'},
                            {'name': 'Reverse reads -', 'ylab': 'Percent converted %'}]
        }        
        
        self.sections.append({
            'name': 'Non T>C mutations over read positions',
            'anchor': 'slamdunk_nontcperreadpos',
            'content': '<p>This plot shows the distribution of non T>C mutations across read positions. <br> \n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperreadpos" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.linegraph.plot([self.nontc_per_readpos_plus, self.nontc_per_readpos_minus], pconfig_nontc) 
        })
        
        self.sections.append({
            'name': 'T>C conversions over read positions',
            'anchor': 'slamdunk_tcperreadpos',
            'content': '<p>This plot shows the distribution of T>C conversions across read positions. <br> \n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperreadpos" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.linegraph.plot([self.tc_per_readpos_plus, self.tc_per_readpos_minus], pconfig_tc) 
        })
        
    def slamdunkTcPerUTRPosPlot (self):
        """ Generate the tc per UTR pos plots """
        
        pconfig_nontc = {
            'id': 'slamdunk_slamdunk_nontcperutrpos_plot',
            'title': 'Non-T>C mutations over UTR ends',
            'ylab': 'Percent mutated %',
            'xlab': 'Position in UTR from 3\' end',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>Pos {point.x}</b>: {point.y:.2f} %',
            'data_labels': [{'name': 'UTRs on plus strand', 'ylab': 'Percent mutated %'},
                            {'name': 'UTRs on minus strand', 'ylab': 'Percent mutated %'}]
        }
        
        pconfig_tc = {
            'id': 'slamdunk_slamdunk_tcperutrpos_plot',
            'title': 'T>C conversions over UTR ends',
            'ylab': 'Percent converted %',
            'xlab': 'Position in UTR from 3\' end',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>Pos {point.x}</b>: {point.y:.2f} %',
            'data_labels': [{'name': 'UTRs on plus strand', 'ylab': 'Percent converted %'},
                            {'name': 'UTRs on minus strand', 'ylab': 'Percent converted %'}]
        }        
        
        self.sections.append({
            'name': 'Non T>C mutations over UTR positions',
            'anchor': 'slamdunk_nontcperutrpos',
            'content': '<p>This plot shows the distribution of non T>C mutations across UTR positions for the last 200 bp from the 3\' UTR end. <br> \n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperutrpos" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.linegraph.plot([self.nontc_per_utrpos_plus, self.nontc_per_utrpos_minus], pconfig_nontc) 
        })
        
        self.sections.append({
            'name': 'T>C conversions over UTR positions',
            'anchor': 'tcperutrpos',
            'content': '<p>This plot shows the distribution of T>C conversions across UTR positions for the last 200 bp from the 3\' UTR end . <br> \n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperutrpos" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.linegraph.plot([self.tc_per_utrpos_plus, self.tc_per_utrpos_minus], pconfig_tc) 
        })