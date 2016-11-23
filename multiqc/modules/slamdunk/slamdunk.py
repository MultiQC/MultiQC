#!/usr/bin/env python

""" MultiQC module to parse output from Slamdunk """

from __future__ import print_function
import logging
import re
from distutils.version import StrictVersion
from collections import OrderedDict

from multiqc import config, BaseMultiqcModule, plots
from curses.ascii import islower

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Slamdunk module class, parses slamdunk logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Slamdunk', anchor='slamdunk',
        href='https://github.com/t-neumann/slamdunk',
        info="is a tool to analyze SLAMSeq data. "\
         "More info to come.")

        self.slamdunk_data = dict()
        self.rates_data_plus = dict()
        self.rates_data_minus = dict()
        self.nontc_per_readpos_plus = dict()
        self.nontc_per_readpos_minus = dict()
        self.tc_per_readpos_plus = dict()
        self.tc_per_readpos_minus = dict()

#         # Find and load any Cutadapt reports
#         self.cutadapt_data = dict()
#         self.cutadapt_length_counts = dict()
#         self.cutadapt_length_exp = dict()
#         self.cutadapt_length_obsexp = dict()
        
        
        for f in self.find_log_files(config.sp['slamdunk']['utrrates']):
            #self.parse_cutadapt_logs(f)
            #log.info(f['f'])
            log.info(f['s_name'])
            log.info(f['root'])
            log.info(f['fn'])
            
        for f in self.find_log_files(config.sp['slamdunk']['summary'], filehandles = True):
            self.parseSummary(f)
            
        for f in self.find_log_files(config.sp['slamdunk']['rates'], filehandles = True):
            self.parseSlamdunkRates(f)
            
        for f in self.find_log_files(config.sp['slamdunk']['tcperreadpos'], filehandles = True):
            self.parseSlamdunkTCPerReadpos(f)

        if len(self.slamdunk_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.slamdunk_data)))
                
        # Start the sections
        self.sections = list()

        # Basic Stats Table
        self.SlamdunkGeneralStatsTable()

        # Rates plot
        self.slamdunkOverallRatesPlot()
        
        # TC per read position plot
        self.slamdunkTcPerReadPosPlot()

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
             
        #log.info(str(baseDict))
        #log.info(str(divisor))
        
        for fromBase in baseDict:
            for toBase in baseDict[fromBase]:
                if(toBase.islower()) :
                    baseDict[fromBase][toBase] = baseDict[fromBase][toBase] / float(divisor[fromBase.lower()]) * 100
                else:
                    baseDict[fromBase][toBase] = baseDict[fromBase][toBase] / float(divisor[fromBase]) * 100
        
        #log.info(str(baseDict))
        
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
    
    def parseSummary(self, f):
        
        # Skip comment line #
        next(f['f'])
        # Skip header line "FileName..."
        next(f['f'])
            
        for line in f['f']:
            fields = line.rstrip().split("\t")
            self.slamdunk_data[self.clean_s_name(fields[0],"")] = dict()
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['sequenced'] = int(fields[4])
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['mapped'] = int(fields[5])
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['deduplicated'] = int(fields[6])
            self.slamdunk_data[self.clean_s_name(fields[0],"")]['filtered'] = int(fields[7])
        self.add_data_source(f)

    def SlamdunkGeneralStatsTable(self):
        """ Take the parsed summary stats from Slamdunk and add it to the
        basic stats table at the top of the report """


        headers = OrderedDict()
        headers['sequenced'] = {
            'title': 'Sequenced',
            'description': '# sequenced reads',
            'shared_key': 'reads',
            'min': 0,
            'format': '{:.f}'
        }
        headers['mapped'] = {
            'title': 'Mapped',
            'description': '# mapped reads',
            'shared_key': 'reads',
            'min': 0,
            'format': '{:.f}'
        }
        headers['deduplicated'] = {
            'title': 'Deduplicated',
            'description': '# deduplicated reads',
            'shared_key': 'reads',
            'min': 0,
            'format': '{:.f}'
        }
        headers['filtered'] = {
            'title': 'Filtered',
            'description': '# reads after filtering',
            'shared_key': 'reads',
            'min': 0,
            'format': '{:.f}'
        }
        self.general_stats_addcols(self.slamdunk_data, headers)
    

    def slamdunkOverallRatesPlot (self):
        """ Generate the overall rates plot """

        pconfig = {
            'id': 'overallratesplot',
            'title': 'Overall conversion rates in reads',
            'cpswitch' : False,
            'cpswitch_c_active': False,
            'stacking' : 'normal',
            'data_labels': [
                "Plus Strand +",
                "Minus Strand -", 
            ]
        }
        
        cats = [OrderedDict(), OrderedDict()]
        
        cats[0]['T>C'] = {
            'name': 'T>C',
            'color': '#D7301F'
        }
        cats[0]['A>T'] = {
            'name': 'A>T',
            'color': '#C6DBEF'
        }
        cats[0]['A>G'] = {
            'name': 'A>G',
            'color': '#6BAED6'
        }
        cats[0]['A>C'] = {
            'name': 'A>C',
            'color': '#2171B5'
        }
        cats[0]['T>A'] = {
            'name': 'T>A',
            'color': '#C7E9C0'
        }
        cats[0]['T>G'] = {
            'name': 'T>G',
            'color': '#74C476'
        }
        cats[0]['G>A'] = {
            'name': 'G>A',
            'color': '#D9D9D9'
        }
        cats[0]['G>T'] = {
            'name': 'G>T',
            'color': '#969696'
        }
        cats[0]['G>C'] = {
            'name': 'G>C',
            'color': '#525252'
        }
        cats[0]['C>A'] = {
            'name': 'C>A',
            'color': '#DADAEB'
        }
        cats[0]['C>T'] = {
            'name': 'C>T',
            'color': '#9E9AC8'
        }
        cats[0]['C>G'] = {
            'name': 'C>G',
            'color': '#6A51A3'
        }
        
        cats[1]['A>G'] = {
            'name': 'A>G',
            'color': '#D7301F'
        }
        cats[1]['A>T'] = {
            'name': 'A>T',
            'color': '#C6DBEF'
        }
        cats[1]['A>C'] = {
            'name': 'A>C',
            'color': '#2171B5'
        }
        cats[1]['T>A'] = {
            'name': 'T>A',
            'color': '#C7E9C0'
        }
        cats[1]['T>G'] = {
            'name': 'T>G',
            'color': '#74C476'
        }
        cats[1]['T>C'] = {
            'name': 'T>C',
            'color': '#238B45'
        }
        cats[1]['G>A'] = {
            'name': 'G>A',
            'color': '#D9D9D9'
        }
        cats[1]['G>T'] = {
            'name': 'G>T',
            'color': '#969696'
        }
        cats[1]['G>C'] = {
            'name': 'G>C',
            'color': '#525252'
        }
        cats[1]['C>A'] = {
            'name': 'C>A',
            'color': '#DADAEB'
        }
        cats[1]['C>T'] = {
            'name': 'C>T',
            'color': '#9E9AC8'
        }
        cats[1]['C>G'] = {
            'name': 'C>G',
            'color': '#6A51A3'
        }
        
        
#         for cat in self.rates_data_plus[self.rates_data_plus.keys()[0]]:
#             cats[0][cat] = dict()
#             cats[0][cat]['name'] = cat
#             cats[0][cat]['color'] = cat
#             
#         for cat in self.rates_data_plus[self.rates_data_minus.keys()[0]]:
#             cats[1][cat] = dict()
#             cats[1][cat]['name'] = cat
        
        self.sections.append({
            'name': 'Conversion rates',
            'anchor': 'overall_rates',
            'content': '<p>This plot shows the individual conversion rates over all reads. \n\
                        It shows these conversion rates strand-specific: This means for a properly labelled \n\
                        sample you would see a T>C excess on \n\
                        the plus-strand and an A>G excess on the minus strand. \n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#rates" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.bargraph.plot([self.rates_data_plus,self.rates_data_minus], cats, pconfig)
        })
        
    def slamdunkTcPerReadPosPlot (self):
        """ Generate the tc per read pos plots """
        
        pconfig_nontc = {
            'id': 'nontcperreadpos_plot',
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
            'id': 'tcperreadpos_plot',
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
            'anchor': 'nontcperreadpos',
            'content': '<p>This plot shows the distribution of non T>C mutations across read positions. \n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperreadpos" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.linegraph.plot([self.nontc_per_readpos_plus, self.nontc_per_readpos_minus], pconfig_nontc) 
        })
        
        self.sections.append({
            'name': 'T>C conversions over read positions',
            'anchor': 'tcperreadpos',
            'content': '<p>This plot shows the distribution of T>C conversions across read positions. \n\
                        See the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperreadpos" target="_blank">slamdunk documentation</a> \n\
                        for more information on how these numbers are generated.</p>' +  
                        plots.linegraph.plot([self.tc_per_readpos_plus, self.tc_per_readpos_minus], pconfig_tc) 
        })