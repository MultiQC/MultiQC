#!/usr/bin/env python

""" MultiQC module to parse output from dragen_markdups """

from __future__ import print_function
from collections import OrderedDict
import logging
from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
	"""dragen_markdups module class, parses stderr logs."""

	def __init__(self):
		# Initialise the parent object
		super(MultiqcModule, self).__init__(name='dragen_markdups', anchor='dragen_markdups',
		info="display the % of duplicate reads from dragen")
		
		# dictionary to hold all data for each sample
		self.dragen_markdups = dict()

		# for each file ending in mappingmetrics.csv
		for f in self.find_log_files('dragen_markdups/mappingmetrics'):
			s_name=f['s_name'].replace(".mapping_metrics",'')
			#take the file name, replace the extension and send this to clean sample name function to get the sample name
			s_name = self.clean_s_name(s_name,f['root'])
			# pass the file and samplename to function which parses file
			parsed_data = self.parse_mapping_metrics(f,s_name)
			# if a result was returned
			if parsed_data is not None:
				if s_name in self.dragen_markdups:
					# write this to log
					log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
				# add the sample as a key to the dictionary and the dictionary of values as the value
				self.dragen_markdups[s_name] = parsed_data[s_name]
				# add data source to multiqc_sources.txt 
				self.add_data_source(f)
		
		# Filter to strip out ignored sample names as per config.yaml
		self.dragen_markdups = self.ignore_samples(self.dragen_markdups)

		if len(self.dragen_markdups) == 0:
			raise UserWarning

		# print number of mapping_metrics files found and parsed
		log.info("Found {} reports".format(len(self.dragen_markdups)))

		# Write parsed report data to a file	
		self.write_data_file(self.dragen_markdups, 'dragen_markdups')

		# add to General Stats Table
		self.dragen_markdups_general_stats_table()


	def parse_mapping_metrics(self, f, s_name):
		""" Go through mapping_metrics.csv file and create a dictionary with the sample name as a key, """
		#create a dictionary to populate from this sample's file
		parsed_data = dict()

		# for each line in the file
		for l in f['f'].splitlines():
			# perform a pattern match to find the line we want 
			if "Number of duplicate reads" in l and "PER RG" not in l:
				# the percentage of duplicate reads is the very last item in the line, seperated by a comma
				# convert this to a float and add to dictionary
				parsed_data[s_name] = {"percentage_duplicates": float(l.split(",")[-1])}
		# else return the dictionary
		return parsed_data

	def dragen_markdups_general_stats_table(self):
		""" Take the percentage_duplicates column and add it to the general stats table at the top of the report"""

		# create a dictionary to hold the columns to add to the general stats table
		headers = OrderedDict()
		# use the column name as per above
		headers['percentage_duplicates'] = {
			'title': '% Dups',
			'description': 'MarkDuplicates - Percent Duplication',
			'format': '{:,.1f}',
			'max': 100,
			'min': 0,
			'suffix': '%',
			'scale': 'OrRd',
			}

		# pass the data dictionary and header dictionary to function to add to table.
		self.general_stats_addcols(self.dragen_markdups, headers)

		