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

		# Adding the bar plots
		self.add_section(
			name = 'Summary',
			anchor = 'seqyclean-summary',
			description = 'This plot shows the number of reads that were kept and discarded',
			plot = self.seqyclean_summary())

		self.add_section(
			name = 'Annotations',
			anchor = 'seqyclean-annotation',
			description = 'This plot shows how reads were annotated',
			plot = self.seqyclean_analysis())

		self.add_section(
			name = 'Discarded',
			anchor = 'seqyclean-discarded',
			description = 'This plot shows the breakdown for discarded reads',
			plot = self.seqyclean_discarded())

		# Write the results to a file
		self.write_data_file(self.seqyclean_data, 'seqyclean')

		# Adding to the general statistics table
		self.seqyclean_general_stats_table()

	def seqyclean_summary(self):
		config = {
			'id': 'seqyclean-1',
			'title': 'SeqyClean: Summary',
			'ylab': 'Number of Reads'
			}
		keys =['PairsKept', # for paired end
			'PairsDiscarded',
			'PE1DiscardedTotal',
			'PE2DiscardedTotal',
			'SEReadsKept', # single end
			'SEDiscardedTotal',
			'ReadsKept', # 454
			'DiscardedTotal']
		return bargraph.plot(self.seqyclean_data, keys, config)

	def seqyclean_analysis(self):
		config = {
			'id': 'seqyclean-2',
			'title': 'SeqyClean: Read Annotations',
			'ylab': 'Number of Reads'
			}
		keys =['PE1TruSeqAdap_found', # for paired end
			'PE1ReadsWVector_found',
			'PE1ReadsWContam_found',
			'PE2TruSeqAdap_found',
			'PE2ReadsWVector_found',
			'PE2ReadsWContam_found',
			'SETruSeqAdap_found', # for single end
			'SEReadsWVector_found',
			'SEReadsWContam_found',
			'left_mid_tags_found', # 454 data
			'right_mid_tags_found',
			'ReadsWithVector_found',
			'ReadsWithContam_found']
		return bargraph.plot(self.seqyclean_data, keys, config)

	def seqyclean_discarded(self):
		config = {
			'id': 'seqyclean-3',
			'title': 'SeqyClean: Discarded Reads',
			'ylab': 'Number of Reads'
			}
		keys =['SEDiscByContam', # single end
			'SEDiscByLength',
			'PE1DiscByContam', # paired end
			'PE1DiscByLength',
			'PE2DiscByContam',
			'PE2DiscByLength',
			'DiscByContam', # 454 data
			'DiscByLength']
		return bargraph.plot(self.seqyclean_data, keys, config)

	def seqyclean_general_stats_table(self):
		headers = OrderedDict()
		# Paired and single end
		headers['Perc_Kept'] = {
			'title': 'Percentage Kept',
			'description': 'The percentage of reads remaining after cleaning',
			'scale': 'YlGn',
			'suffix': '%',
			'max': 100,
			'min': 0
		}
		# 454
		headers['PercentageKept'] = {
			'title': 'Percentage Kept',
			'description': 'The percentage of reads remaining after cleaning',
			'scale': 'YlGn',
			'suffix': '%',
			'max': 100,
			'min': 0
		}
		self.general_stats_addcols(self.seqyclean_data, headers)
