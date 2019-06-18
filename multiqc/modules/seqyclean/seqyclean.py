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
		info="Analysizes the seqyclean files in the pipeline.")

		# Parse logs
		self.seqyclean = dict()
		data = {}
		for f in self.find_log_files('seqyclean'):
			RowTwo = f['f'].split("\n")[1] #split the file into two rows
			PairsKept = RowTwo.split("\t")[57] # get the number of pairs kept (from the second row)
			PairsDiscarded = RowTwo.split("\t")[61] #get num of pairs discarded
			data.update({
				f['s_name']: {
					'PairsKept': PairsKept,
					'PairsDiscarded': PairsDiscarded
				}
			})
		# Write the results to a file
		self.write_data_file(data, 'seqyclean_results')
		# Create files
		headers = OrderedDict()
		headers['PairsKept'] = {
			'title': 'PairsKept',
			'description': 'Num of pairs kept from seqyclean analysis',
			'scale': 'YlGn',
			'format': '{} bp'
		}
		headers['PairsDiscarded'] = {
			'title': 'PairsDiscarded',
			'description': 'Num of pairs discarded from seqyclean analysis',
			'scale': 'OrRd',
			'format': '{} bp'
		}
		self.general_stats_addcols(data, headers)
		
		self.add_section (
			name = 'Seqyclean Results',
			anchor = 'seqyclean-first',
			description = 'This shows the results of the seqyclean process',
			helptext = "If you're not sure _how_ to interpret the data, we can help!",
			plot = bargraph.plot(data)
		)

	"""	
	#Used this originally but now I don't have to
	def parse_seqyclean(self, f):
		#l = f['s_name'].split("-",3)
		RowTwo = f['f'].split("\n")[1] #split the file into two rows
		PairsKept = RowTwo.split("\t")[57] # get the number of pairs kept (from the second row)
		PairsDiscarded = RowTwo.split("\t")[61] #get num of pairs discarded
		parsed_data = {
			f['s_name'] : (PairsKept, PairsDiscarded)
		}
	
		return parsed_data
		"""
	#print( f )       # File contents
	#print (f['f]) This just prints the contents of the file, and not the name, root, etc..
	#print( f['s_name'] )  # Sample name (from cleaned filename)
	#print( f['fn'] )      # Filename
	#print( f['root'] ) 
	
'''#!/usr/bin/env python3

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
		info="Analysizes the seqyclean files in the pipeline.")

		# Parse logs
		self.seqyclean = dict()
		data = {}
		for f in self.find_log_files('seqyclean'):
			data.update(self.parse_seqyclean(f))
		#self.write_data_file(data, 'seqyclean')
		headers = OrderedDict()
		headers['PairsKept'] = {
			'title': 'PairsKept',
			'description': 'Num of pairs kept from seqyclean analysis',
			'scale': 'PuBuGn'
		}
		headers['PairsDiscarded'] = {
			'title': 'PairsDiscarded',
			'description': 'Num of pairs discarded from seqyclean analysis',
			'scale': 'PiYG'
		}

		#self.general_stats_addcols(data, headers)
		
	def parse_seqyclean(self, f):
		#l = f['s_name'].split("-",3)
		RowTwo = f['f'].split("\n")[1] #split the file into two rows
		PairsKept = RowTwo.split("\t")[57] # get the number of pairs kept (from the second row)
		PairsDiscarded = RowTwo.split("\t")[61] #get num of pairs discarded
		parsed_data = {
			f['s_name'] : (PairsKept, PairsDiscarded)
		}
	
		return parsed_data
	
		
	#print( f )       # File contents
	#print (f['f]) This just prints the contents of the file, and not the name, root, etc..
	#print( f['s_name'] )  # Sample name (from cleaned filename)
	#print( f['fn'] )      # Filename
	#print( f['root'] ) 
	'''