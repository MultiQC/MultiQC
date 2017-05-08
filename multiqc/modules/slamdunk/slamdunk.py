#!/usr/bin/env python

""" MultiQC module to parse output from Slamdunk """

from __future__ import print_function
import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table, bargraph, linegraph, scatter

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

        num_reports = 0

        self.plot_cols = [
            '#fdbf6f','#2171B5','#6BAED6','#C6DBEF','#74C476','#C7E9C0',
            '#D9D9D9','#969696','#525252','#DADAEB','#9E9AC8','#6A51A3'
        ]

        # Summary Reports
        self.slamdunk_data = dict()
        for f in self.find_log_files('slamdunk/summary', filehandles = True):
            self.parseSummary(f)
        self.slamdunk_data = self.ignore_samples(self.slamdunk_data)
        if len(self.slamdunk_data) > 0:
            self.slamdunkGeneralStatsTable()
            self.slamdunkFilterStatsTable()
            log.debug("Found {} summary reports".format(len(self.slamdunk_data)))
            num_reports = max(num_reports, len(self.slamdunk_data))

        # PCA Plots
        self.PCA_data = dict()
        for f in self.find_log_files('slamdunk/PCA', filehandles = True):
            self.parsePCA(f)
        self.PCA_data = self.ignore_samples(self.PCA_data)
        if len(self.PCA_data) > 0:
            self.slamdunkPCAPlot()
            log.debug("Found {} PCA plots".format(len(self.PCA_data)))
            num_reports = max(num_reports, len(self.PCA_data))

        # UTR Rate reports
        self.utrates_data = dict()
        for f in self.find_log_files('slamdunk/utrrates', filehandles = True):
            self.parseUtrRates(f)
        self.utrates_data = self.ignore_samples(self.utrates_data)
        if len(self.utrates_data) > 0:
            self.write_data_file(self.utrates_data, 'multiqc_slamdunk_utrrates')
            self.slamdunkUtrRatesPlot()
            log.debug("Found {} UTR rate reports".format(len(self.utrates_data)))
            num_reports = max(num_reports, len(self.utrates_data))

        # Read rate reports
        self.rates_data_plus = dict()
        self.rates_data_minus = dict()
        for f in self.find_log_files('slamdunk/rates', filehandles = True):
            self.parseSlamdunkRates(f)
        self.rates_data_plus = self.ignore_samples(self.rates_data_plus)
        self.rates_data_minus = self.ignore_samples(self.rates_data_minus)
        if len(self.rates_data_plus) > 0:
            self.write_data_file(self.rates_data_plus, 'multiqc_slamdunk_readrates_plus')
            self.write_data_file(self.rates_data_minus, 'multiqc_slamdunk_readrates_minus')
            self.slamdunkOverallRatesPlot()
            log.debug("Found {} read rate reports".format(len(self.rates_data_plus)))
            num_reports = max(num_reports, len(self.rates_data_plus))

        # TCP error read rate
        self.nontc_per_readpos_plus = dict()
        self.nontc_per_readpos_minus = dict()
        self.tc_per_readpos_plus = dict()
        self.tc_per_readpos_minus = dict()
        for f in self.find_log_files('slamdunk/tcperreadpos', filehandles = True):
            self.parseSlamdunkTCPerReadpos(f)
        self.nontc_per_readpos_plus = self.ignore_samples(self.nontc_per_readpos_plus)
        self.nontc_per_readpos_minus = self.ignore_samples(self.nontc_per_readpos_minus)
        self.tc_per_readpos_plus = self.ignore_samples(self.tc_per_readpos_plus)
        self.tc_per_readpos_minus = self.ignore_samples(self.tc_per_readpos_minus)
        if len(self.tc_per_readpos_plus) > 0:
            self.write_data_file(self.tc_per_readpos_plus, 'multiqc_slamdunk_tcperreadpos_plus')
            self.write_data_file(self.nontc_per_readpos_plus, 'multiqc_slamdunk_nontcperreadpos_plus')
            self.write_data_file(self.tc_per_readpos_minus, 'multiqc_slamdunk_tcperreadpos_minus')
            self.write_data_file(self.nontc_per_readpos_minus, 'multiqc_slamdunk_nontcperreadpos_minus')
            self.slamdunkTcPerReadPosPlot()
            log.debug("Found {} TCP error read rate reports".format(len(self.tc_per_readpos_plus)))
            num_reports = max(num_reports, len(self.tc_per_readpos_plus))

        # Non-TCP error read rate
        self.nontc_per_utrpos_plus = dict()
        self.nontc_per_utrpos_minus = dict()
        self.tc_per_utrpos_plus = dict()
        self.tc_per_utrpos_minus = dict()
        for f in self.find_log_files('slamdunk/tcperutrpos', filehandles = True):
            self.parseSlamdunkTCPerUtrpos(f)
        self.nontc_per_utrpos_plus = self.ignore_samples(self.nontc_per_utrpos_plus)
        self.nontc_per_utrpos_minus = self.ignore_samples(self.nontc_per_utrpos_minus)
        self.tc_per_utrpos_plus = self.ignore_samples(self.tc_per_utrpos_plus)
        self.tc_per_utrpos_minus = self.ignore_samples(self.tc_per_utrpos_minus)
        if len(self.nontc_per_utrpos_plus) > 0:
            self.write_data_file(self.tc_per_utrpos_plus, 'multiqc_slamdunk_tcperutrpos_plus')
            self.write_data_file(self.nontc_per_utrpos_plus, 'multiqc_slamdunk_nontcperutrpos_plus')
            self.write_data_file(self.tc_per_utrpos_minus, 'multiqc_slamdunk_tcperutrpos_minus')
            self.write_data_file(self.nontc_per_utrpos_minus, 'multiqc_slamdunk_nontcperutrpos_minus')
            self.slamdunkTcPerUTRPosPlot()
            log.debug("Found {} non TCP error read rate reports".format(len(self.nontc_per_utrpos_plus)))
            num_reports = max(num_reports, len(self.nontc_per_utrpos_plus))

        if num_reports == 0:
            log.debug("No slamdunk reports found.")
            raise UserWarning
        else:
            log.info("Found {} reports".format(num_reports))


    def parsePCA(self, f):

        # Skip header
        next(f['f'])

        for line in f['f']:
            fields = line.rstrip().split('\t')

            sample = self.clean_s_name(fields[0],f['root'])
            PC1 = fields[1]
            PC2 = fields[2]

            self.PCA_data[sample] = [{'x': float(PC1), 'y': float(PC2)}]


    def parseUtrRates(self, f):

        # Skip comment line #
        next(f['f'])

        # Read median header
        line = next(f['f'])

        if "Conversions=" in line:

            sample = f['s_name']
            self.utrates_data[sample] = OrderedDict()

            conversions = re.sub(".*Conversions=","",line.rstrip()).split(",")

            for conversion in conversions:
                type, value = conversion.split(":")
                self.utrates_data[sample][type] = float(value)

        else:
            log.warning("Malformed UTR rates header. Conversion rates per UTR plot will be affected.")

    def parseSlamdunkRates(self, f):

        sample = f['s_name']

        # Skip comment line #
        next(f['f'])

        bases = next(f['f']).rstrip().split('\t')

        baseDict = {}
        order = {}

        for i in range(1, len(bases)):
            order[i] = bases[i]

        for line in f['f']:
            values = line.rstrip().split('\t')
            base = values[0]
            baseDict[base]= {}

            for i in range(1, len(values)):
                baseDict[base][order[i]] = int(values[i])

        divisor = {}

        for fromBase in baseDict:
            for toBase in baseDict[fromBase]:
                if(toBase.islower()):
                    if not fromBase.lower() in divisor:
                        divisor[fromBase.lower()] = 0
                    divisor[fromBase.lower()] += baseDict[fromBase][toBase]
                else:
                    if not fromBase in divisor:
                        divisor[fromBase] = 0
                    divisor[fromBase] += baseDict[fromBase][toBase]

        for fromBase in baseDict:
            for toBase in baseDict[fromBase]:
                if(toBase.islower()):
                    if divisor[fromBase.lower()] > 0:
                        baseDict[fromBase][toBase] = baseDict[fromBase][toBase] / float(divisor[fromBase.lower()]) * 100
                    else :
                        baseDict[fromBase][toBase] = 0.0
                else:
                    if divisor[fromBase] > 0:
                        baseDict[fromBase][toBase] = baseDict[fromBase][toBase] / float(divisor[fromBase]) * 100
                    else :
                        baseDict[fromBase][toBase] = 0.0

        self.rates_data_plus[sample] = {}
        self.rates_data_minus[sample] = {}

        for fromBase in baseDict:
            for toBase in baseDict[fromBase]:
                if fromBase != "N" and toBase.upper() != "N" and fromBase != toBase.upper():
                    if(toBase.islower()):
                        self.rates_data_minus[sample][fromBase + ">" + toBase.upper()] = baseDict[fromBase][toBase]
                    else:
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
            if int(values[4]) > 0:
                self.nontc_per_readpos_plus[sample][pos] = float(int(values[0])) / int(values[4]) * 100
                self.tc_per_readpos_plus[sample][pos] = float(int(values[2])) / int(values[4]) * 100
            else:
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
            if int(values[4]) > 0:
                self.nontc_per_utrpos_plus[sample][pos] = float(int(values[0])) / int(values[4]) * 100
                self.tc_per_utrpos_plus[sample][pos] = float(int(values[2])) / int(values[4]) * 100
            else:
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
            s_name = self.clean_s_name(fields[0],f['root'])
            self.slamdunk_data[s_name] = dict()
            self.slamdunk_data[s_name]['sequenced'] = int(fields[4])
            self.slamdunk_data[s_name]['mapped'] = int(fields[5])
            #self.slamdunk_data[s_name]['deduplicated'] = int(fields[6])
            self.slamdunk_data[s_name]['mqfiltered'] = int(fields[7])
            self.slamdunk_data[s_name]['idfiltered'] = int(fields[8])
            self.slamdunk_data[s_name]['nmfiltered'] = int(fields[9])
            self.slamdunk_data[s_name]['multimapper'] = int(fields[10])
            self.slamdunk_data[s_name]['retained'] = int(fields[11])

            # Additional Count Column found in Table
            if columnCount == 14:
                self.slamdunk_data[s_name]['counted'] = int(fields[12])

        self.add_data_source(f)

    def slamdunkGeneralStatsTable(self):
        """ Take the parsed summary stats from Slamdunk and add it to the
        basic stats table at the top of the report """


        headers = OrderedDict()
        headers['counted'] = {
            'title': '{} Counted'.format(config.read_count_prefix),
            'description': '# reads counted within 3\'UTRs ({})'.format(config.read_count_desc),
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * config.read_count_multiplier,
        }
        headers['retained'] = {
            'title': '{} Retained'.format(config.read_count_prefix),
            'description': '# retained reads after filtering ({})'.format(config.read_count_desc),
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * config.read_count_multiplier,
        }
        headers['mapped'] = {
            'title': '{} Mapped'.format(config.read_count_prefix),
            'description': '# mapped reads ({})'.format(config.read_count_desc),
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * config.read_count_multiplier,
        }
        headers['sequenced'] = {
            'title': '{} Sequenced'.format(config.read_count_prefix),
            'description': '# sequenced reads ({})'.format(config.read_count_desc),
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * config.read_count_multiplier,
        }

        self.general_stats_addcols(self.slamdunk_data, headers)

    def slamdunkFilterStatsTable(self):
        """ Take the parsed filter stats from Slamdunk and add it to a separate table """

        headers = OrderedDict()
        headers['mapped'] = {
            'namespace': 'Slamdunk',
            'title': '{} Mapped'.format(config.read_count_prefix),
            'description': '# mapped reads ({})'.format(config.read_count_desc),
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:,.2f}',
            'suffix': config.read_count_prefix,
            'scale': 'YlGn',
            'modify': lambda x: float(x) * config.read_count_multiplier,
        }
        headers['multimapper'] = {
            'namespace': 'Slamdunk',
            'title': '{} Multimap-Filtered'.format(config.read_count_prefix),
            'description': '# multimap-filtered reads ({})'.format(config.read_count_desc),
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:,.2f}',
            'suffix': config.read_count_prefix,
            'scale': 'OrRd',
            'modify': lambda x: float(x) * config.read_count_multiplier,
        }
        headers['nmfiltered'] = {
            'namespace': 'Slamdunk',
            'title': '{} NM-Filtered'.format(config.read_count_prefix),
            'description': '# NM-filtered reads ({})'.format(config.read_count_desc),
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:,.2f}',
            'suffix': config.read_count_prefix,
            'scale': 'OrRd',
            'modify': lambda x: float(x) * config.read_count_multiplier,
        }
        headers['idfiltered'] = {
            'namespace': 'Slamdunk',
            'title': '{} Identity-Filtered'.format(config.read_count_prefix),
            'description': '# identity-filtered reads ({})'.format(config.read_count_desc),
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:,.2f}',
            'suffix': config.read_count_prefix,
            'scale': 'OrRd',
            'modify': lambda x: float(x) * config.read_count_multiplier,
        }
        headers['mqfiltered'] = {
            'namespace': 'Slamdunk',
            'title': '{} MQ-Filtered'.format(config.read_count_prefix),
            'description': '# MQ-filtered reads ({})'.format(config.read_count_desc),
            'shared_key': 'read_count',
            'min': 0,
            'format': '{:,.2f}',
            'suffix': config.read_count_prefix,
            'scale': 'OrRd',
            'modify': lambda x: float(x) * config.read_count_multiplier,
        }
        pconfig = {
            'id': 'slamdunk_filtering_table',
            'min': 0,
        }

        self.add_section (
            name = 'Filter statistics',
            anchor = 'slamdunk_filtering',
            description = 'This table shows the number of reads filtered with each filter criterion during filtering phase of slamdunk.',
            plot = table.plot(self.slamdunk_data, headers, pconfig)
        )

    def slamdunkOverallRatesPlot (self):
        """ Generate the overall rates plot """

        pconfig = {
            'id': 'overallratesplot',
            'title': 'Overall conversion rates in reads',
            'cpswitch': False,
            'cpswitch_c_active': False,
            'stacking': 'normal',
            'tt_decimals': 2,
            'tt_suffix': '%',
            'tt_percentages': False,
            'hide_zero_cats': False,
            'data_labels': [
                "Plus Strand +",
                "Minus Strand -",
            ]
        }

        cats = [OrderedDict(), OrderedDict()]
        keys = [
            ['T>C', 'A>T', 'A>G', 'A>C', 'T>A', 'T>G', 'G>A', 'G>T', 'G>C', 'C>A', 'C>T', 'C>G'],
            ['A>G','A>T','A>C','T>A','T>G','T>C','G>A','G>T','G>C','C>A','C>T','C>G']
        ]
        for i, k in enumerate(keys):
            for j, v in enumerate(k):
                cats[i][v] = { 'color': self.plot_cols[j] }

        self.add_section (
            name = 'Conversion rates per read',
            anchor = 'slamdunk_overall_rates',
            description = """This plot shows the individual conversion rates over all reads.
                        It shows these conversion rates strand-specific: This means for a properly labelled
                        sample you would see a T&gt;C excess on the plus-strand and an A&gt;G excess on the minus strand
                        (see the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#rates" target="_blank">slamdunk docs</a>).""",
            plot = bargraph.plot([self.rates_data_plus,self.rates_data_minus], cats, pconfig)
        )


    def slamdunkUtrRatesPlot (self):
        """ Generate the UTR rates plot """

        cats = OrderedDict()
        keys = ['T>C', 'A>T', 'A>G', 'A>C', 'T>A', 'T>G', 'G>A', 'G>T', 'G>C', 'C>A', 'C>T', 'C>G']
        for i, v in enumerate(keys):
            cats[v] = { 'color': self.plot_cols[i] }

        pconfig = {
            'id': 'slamdunk_utrratesplot',
            'title': 'Overall conversion rates per UTR',
            'cpswitch': False,
            'cpswitch_c_active': False,
            'stacking': 'normal',
            'tt_decimals': 2,
            'tt_suffix': '%',
            'tt_percentages': False,
            'hide_zero_cats': False
        }

        self.add_section (
            name = 'Conversion rates per UTR',
            anchor = 'slamdunk_utr_rates',
            description = """This plot shows the individual conversion rates for all UTRs
                        (see the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#utrrates" target="_blank">slamdunk docs</a>).""",
            plot = bargraph.plot(self.utrates_data, cats, pconfig)
        )

    def slamdunkPCAPlot (self):
        """ Generate the PCA plots """

        pconfig = {
            'id': 'slamdunk_pca',
            'title': 'Slamdunk PCA',
            'xlab': 'PC1',
            'ylab': 'PC2',
            'tt_label': 'PC1 {point.x:.2f}: PC2 {point.y:.2f}'
        }

        self.add_section (
            name = 'PCA (T&gt;C based)',
            anchor = 'slamdunk_PCA',
            description = """This plot shows the principal components of samples based
                        on the distribution of reads with T&gt;C conversions within UTRs
                        (see the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#summary" target="_blank">slamdunk docs</a>).""",
            plot = scatter.plot(self.PCA_data, pconfig)
        )

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

        self.add_section (
            name = 'Non T&gt;C mutations over read positions',
            anchor = 'slamdunk_nontcperreadpos',
            description = """This plot shows the distribution of non T&gt;C mutations across read positions
                        (see the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperreadpos" target="_blank">slamdunk docs</a>).""",
            plot = linegraph.plot([self.nontc_per_readpos_plus, self.nontc_per_readpos_minus], pconfig_nontc)
        )

        self.add_section (
            name = 'T&gt;C conversions over read positions',
            anchor = 'slamdunk_tcperreadpos',
            description = """This plot shows the distribution of T&gt;C conversions across read positions
                        (see the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperreadpos" target="_blank">slamdunk docs</a>).""",
            plot = linegraph.plot([self.tc_per_readpos_plus, self.tc_per_readpos_minus], pconfig_tc)
        )

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

        self.add_section (
            name = 'Non T&gt;C mutations over UTR positions',
            anchor = 'slamdunk_nontcperutrpos',
            description = """This plot shows the distribution of non T&gt;C mutations across UTR positions for the last 200 bp from the 3\' UTR end
                        (see the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperutrpos" target="_blank">slamdunk docs</a>).""",
            plot = linegraph.plot([self.nontc_per_utrpos_plus, self.nontc_per_utrpos_minus], pconfig_nontc)
        )

        self.add_section (
            name = 'T&gt;C conversions over UTR positions',
            anchor = 'tcperutrpos',
            description = """This plot shows the distribution of T&gt;C conversions across UTR positions for the last 200 bp from the 3\' UTR end
                        (see the <a href="http://slamdunk.readthedocs.io/en/latest/Alleyoop.html#tcperutrpos" target="_blank">slamdunk docs</a>).""",
            plot = linegraph.plot([self.tc_per_utrpos_plus, self.tc_per_utrpos_minus], pconfig_tc)
        )
