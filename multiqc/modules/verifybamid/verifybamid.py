#!/usr/bin/env python

""" MultiQC module to parse output from Peddy """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Peddy module class, parses stderr logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Verify BAM ID', anchor='verifybamid',
        href='https://genome.sph.umich.edu/wiki/VerifyBamID',
        info="blah")
        
        # dictionary to hold all data for each sample
        self.verify_data = dict()

        # for each file ending in self.SM
        for f in self.find_log_files('verifybamid/selfsm'):
            # pass the file to function self.parse_selfsm to parse file
            parsed_data = self.parse_selfsm(f)
            # if data was found clean the sample name
            if parsed_data is not None:
                for s_name in parsed_data:
                    # clean the sample name
                    s_name = self.clean_s_name(s_name, f['root'])
                    try:
                        self.verify_data[s_name].update(parsed_data[s_name])
                    except KeyError:
                        self.verify_data[s_name] = parsed_data[s_name]

        # Filter to strip out ignored sample names
        self.verify_data = self.ignore_samples(self.verify_data)

        if len(self.verify_data) == 0:
            raise UserWarning

        # print number of verifyBAMID reports found and parsed
        log.info("Found {} reports".format(len(self.verify_data)))

        # Write parsed report data to a file	
        self.write_data_file(self.verify_data, 'multiqc_verifybamid')

        # Basic Stats Table
        self.peddy_general_stats_table()


def parse_selfsm(self, f):
    """ Go through selfSM file and create a dictionary with the sample name as a key, """
    #create a dictionary to populate from this sample's file
    parsed_data = dict()
    # set a empty variable which denotes if the headers have been read
    headers = None
    # for each line in the file
    for l in f['f'].splitlines():
            # split the line on tab
            s = l.split("\t")
            # if we haven't already read the header line
            if headers is None:
                # remove the leading hash from the first element of the header
                s[0] = s[0].lstrip('#')
                # assign this list to headers variable
                headers = s
            # for all rows after the first
            else:
                # create a dictionary entry with the first column as a key (sample name) and an empty dictionary as a value
                parsed_data[s[0]] = dict()
                # for each item in list of items in the row
                for i, v in enumerate(s):
                # if it's not the first element (if it's not the name)
                    if i != 0:
                        # try and convert the value into a float
                        try:
                        # and add to the dictionary the key as the corrsponding item from the header and the value from the list
                            parsed_data[s[0]][headers[i]] = float(v)
                            #if can't convert to float...
                        except ValueError:
                        # add to the dictionary the key as the corrsponding item from the header and the value from the list
                            parsed_data[s[0]][headers[i]] = v
    # nothing has been parsed return None
    if len(parsed_data) == 0:
        return None
    # else return the dictionary
    return parsed_data

    def peddy_general_stats_table(self):
        """ Take the percentage of contamination from all the parsed *.SELFSM files and add it to the basic stats table at the top of the report """

        # create a dictionary to hold the columns to add to the general stats table
        headers = OrderedDict()

        # add the FREEMIX column. set the title and description
        headers['FREEMIX'] = {
        	'title': 'contamination_prediction',
        	'description': 'VerifyBamID predicted level of contamination',
        }
        # pass the data dictionary and header dictionary to function to add to table.
        self.general_stats_addcols(self.verify_data, headers)

	
