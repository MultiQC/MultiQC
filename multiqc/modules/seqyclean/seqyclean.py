# !/usr/bin/env python

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import linegraph, bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
	def __init__(self):
		# Initialise the parent object
		super(MultiqcModule, self).__init__(name='SeqyClean',
		anchor='seqyclean',
		href='https://github.com/ibest/seqyclean',
		info='SeqyClean is a pre-processing tool for NGS data that filters adapters, vectors, and contaminants while quality trimming.')

		# Parse logs
		self.seqyclean_data = dict()
		for f in self.find_log_files('seqyclean'):
			rows = f['f'].splitlines()
			headers = rows[0].split("\t")
			cols = rows[1].split("\t")

			self.seqyclean_data[f['s_name']] = dict()
			for header, cols in zip(headers, cols):
				self.seqyclean_data[f['s_name']].update({ header : cols })

		if len(self.seqyclean_data) == 0:
			raise UserWarning

		self.seqyclean_data = self.ignore_samples(self.seqyclean_data)

		log.info("Found {} logs".format(len(self.seqyclean_data)))

		# Adding the bar plot
		self.add_section(plot = self.seqyclean_bargraph())

		# Write the results to a file
		self.write_data_file(self.seqyclean_data, 'seqyclean')

		# Adding to the general statistics table
		self.seqyclean_general_stats_table()

	def seqyclean_bargraph(self):
		config = {
			'id': 'seqyclean-1',
			'title': 'SeqyClean: Reads Analysis',
			'ylab': 'Number of Reads'
			}
		keys =['PE1ReadsAn',
			'PE1TruSeqAdap_found',
			'PE1ReadswVector_found',
			'PE1ReadswContam_found',
			'PE1DiscByContam',
			'PE1DiscByLength',
			'PE2ReadsAn',
			'PE2TruSeqAdap_found',
			'PE2ReadswVector_found',
			'PE2ReadswContam_found',
			'PE2DiscByContam',
			'PE2DiscByLength']

		return bargraph.plot(self.seqyclean_data, keys, config)

	def seqyclean_general_stats_table(self):
		headers = OrderedDict()
		headers['PairsKept'] = {
			'title': '{} Pairs Kept'.format(config.read_count_prefix),
			'description': 'The number of read pairs remaining ({})'.format(config.read_count_desc),
			'modify': lambda x: x * config.read_count_multiplier,
			'shared_key': 'read_count',
			'scale': 'YlGn',
			'id': 'pairs_kept',
			'format' : '{:,.0f}',
		}
		headers['PairsDiscarded'] = {
			'title': 'Pairs Discarded',
			'description': 'The number of read pairs discarded ({})'.format(config.read_count_desc),
			'modify': lambda x: x * config.read_count_multiplier,
			'shared_key': 'read_count',
			'scale': 'OrRd',
			'format' : '{:,.0f}',
			'id': 'pairs_discarded'
		}
		self.general_stats_addcols(self.seqyclean_data, headers)
