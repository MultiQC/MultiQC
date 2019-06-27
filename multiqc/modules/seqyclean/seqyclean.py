#!/usr/bin/env python3

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import linegraph, bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)
import traceback

class MultiqcModule(BaseMultiqcModule):
	def __init__(self):
		# Initialise the parent object
		super(MultiqcModule, self).__init__(name='Seqyclean', anchor='seqyclean',
		href="http://www.awesome_bioinfo.com/my_module",
		info="These graphs analyze the seqyclean files created during the pipeline. There are two associated graphs. The first one shows the total number of pairs kept and discarded. The second breaks down the reasons why certain pairs were discarded.")

		# Parse logs
		self.seqyclean = dict()
		self.data_total = {}
		self.data_breakdown = {}
		for f in self.find_log_files('seqyclean'):
			RowTwo = f['f'].split("\n")[1] #split the file into two rows
			#These variables will be used in the total data table 
			PairsKept = RowTwo.split("\t")[57] # get the number of pairs kept (from the second row)
			PairsDiscarded = RowTwo.split("\t")[61] #get num of pairs discarded
			#These are for the breakdown
			PE1ReadsAn = RowTwo.split("\t")[21]
			PE1TruSeqAdap_found = RowTwo.split("\t")[23]
			PE1ReadswVector_found = RowTwo.split("\t")[25]
			PE1ReadswContam_found = RowTwo.split("\t")[27]
			PE1DiscByContam = RowTwo.split("\t")[37]
			PE1DiscByLength = RowTwo.split("\t")[38]
			PE2ReadsAn = RowTwo.split("\t")[39]
			PE2TruSeqAdap_found = RowTwo.split("\t")[41]
			PE2ReadswVector_found = RowTwo.split("\t")[43]
			PE2ReadswContam_found = RowTwo.split("\t")[45]
			PE2DiscByContam = RowTwo.split("\t")[55]
			PE2DiscByLength = RowTwo.split("\t")[56]

			#add the variables above into a dictionary for the total values table
			self.data_total.update({
				f['s_name']: {
					'PairsKept': PairsKept,
					'PairsDiscarded': PairsDiscarded
				}
			})
			#add the above variables into a dictionary for a breakdown table
			self.data_breakdown.update({
				f['s_name']: {
					'PE1ReadsAn':PE1ReadsAn,
					'PE1TruSeqAdap_found':PE1TruSeqAdap_found,
					'PE1ReadswVector_found':PE1ReadswVector_found,
					'PE1ReadswContam_found':PE1ReadswContam_found,
					'PE1DiscByContam':PE1DiscByContam,
					'PE1DiscByLength':PE1DiscByLength,
					'PE2ReadsAn':PE2ReadsAn,
					'PE2TruSeqAdap_found':PE2TruSeqAdap_found,
					'PE2ReadswVector_found':PE2ReadswVector_found,
					'PE2ReadswContam_found':PE2ReadswContam_found,
					'PE2DiscByContam':PE2DiscByContam,
					'PE2DiscByLength':PE2DiscByLength
				}
			})
		
		# Write the results to a file
		self.write_data_file(self.data_total, 'seqyclean_results')
		self.write_data_file(self.data_breakdown, 'seqyclean_results')
		#headers for the "general information" at the top of the Multiqc Analysis
		self.seqyclean_general_stats_table()
		#These functions create the sections in the report summary
		firstkeys = ['PairsKept', 'PairsDiscarded']
		secondkeys =['PE1ReadsAn','PE1TruSeqAdap_found','PE1ReadswVector_found','PE1ReadswContam_found','PE1DiscByContam','PE1DiscByLength','PE2ReadsAn','PE2TruSeqAdap_found','PE2ReadswVector_found','PE2ReadswContam_found','PE2DiscByContam','PE2DiscByLength']
		pconfig = {
            'id': 'seqyclean_assignment_plot',
            'title': 'Seqyclean: Reads',
            'ylab': 'Num of Reads',
        }
		config = {
            'id': 'seqyclean_assignment_plot_2',
            'title': 'Seqyclean: Breakdown Reads',
            'ylab': 'Num of Reads',
        }
		self.add_section (
			name = 'Total Seqyclean Results',
			anchor = 'seqyclean-first',
			description = 'This shows the results of the seqyclean process',
			helptext = "If you're not sure _how_ to interpret the data, we can help!",
			plot = bargraph.plot(self.data_total, firstkeys, pconfig)
		)
		self.add_section (
			name = 'Seqyclean Results Breakdown',
			anchor = 'seqyclean-second',
			description = 'This shows the breakdown results of the seqyclean process',
			plot = bargraph.plot(self.data_breakdown, secondkeys, config)
		)
				
	def seqyclean_general_stats_table(self):
		""" Take the parsed stats from the sargasso report and add them to the
		basic stats table at the top of the report """
		headers = OrderedDict()
		headers['PairsKept'] = {
			'title': 'PairsKept: Seqyclean pairs that were preserved',
			'description': 'Num of pairs kept from seqyclean analysis',
			'scale': 'YlGn',
			'format': '{} bp',
			'ylab': '# Reads',
			'id': 'pairs_kept'
		}
		headers['PairsDiscarded'] = {
			'title': 'PairsDiscarded: Pairs that were discarded',
			'description': 'Num of pairs discarded from seqyclean analysis',
			'scale': 'OrRd',
			'format': '{} bp',
			'ylab': '# Reads',
			'id': 'pairs_discarded'
		}
		self.general_stats_addcols(self.data_total, headers)