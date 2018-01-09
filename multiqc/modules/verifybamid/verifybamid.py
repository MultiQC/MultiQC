#!/usr/bin/env python

""" MultiQC module to parse output from VerifyBAMID """

from __future__ import print_function
from collections import OrderedDict
import logging
from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
     module class, parses stderr logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='verifyBAMID', anchor='verifybamid',
        href='https://genome.sph.umich.edu/wiki/VerifyBamID',
        info="detect sample contamination and/or sample swap")
        
        # dictionary to hold all data for each sample
        self.verifybamid_data = dict()

        # for each file ending in self.SM
        for f in self.find_log_files('verifybamid/selfsm'):
            # add data source to multiqc_sources.txt 
            self.add_data_source(f)
            # pass the file to function self.parse_selfsm to parse file
            parsed_data = self.parse_selfsm(f)
            # if a result was returned
            if parsed_data is not None:
                # for each sample extracted from the file
                for s_name in parsed_data:
                    # add the sample as a key to the verifybamid_data dictionary and the dictionary of values as the value
                    self.verifybamid_data[s_name] = parsed_data[s_name]
        
        # Filter to strip out ignored sample names
        self.verifybamid_data = self.ignore_samples(self.verifybamid_data)

        if len(self.verifybamid_data) == 0:
            raise UserWarning

        # print number of verifyBAMID reports found and parsed
        log.info("Found {} reports".format(len(self.verifybamid_data)))

        # Write parsed report data to a file	
        self.write_data_file(self.verifybamid_data, 'multiqc_verifybamid')

        # add to General Stats Table
        self.verifybamid_general_stats_table()

        # add section with the values from the verify BAM ID output
        self.verifybamid_table()


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
                    s_name=s[0]
                    s_name = self.clean_s_name(s_name, f['root'])
                    # create a dictionary entry with the first column as a key (sample name) and an empty dictionary as a value
                    parsed_data[s_name] = {}
                    # for each item in list of items in the row
                    for i, v in enumerate(s):
                    # if it's not the first element (if it's not the name)
                        if i != 0:
                            # try and convert the value into a float
                            try:
                            # and add to the dictionary the key as the corrsponding item from the header and the value from the list
                                parsed_data[s_name][headers[i]] = float(v)
                                #if can't convert to float...
                            except ValueError:
                            # add to the dictionary the key as the corrsponding item from the header and the value from the list
                                parsed_data[s_name][headers[i]] = v
    
        # nothing has been parsed return None
        if len(parsed_data) == 0:
            return None
        # else return the dictionary
        return parsed_data

    def verifybamid_general_stats_table(self):
        """ Take the percentage of contamination from all the parsed *.SELFSM files and add it to the basic stats table at the top of the report """

        # create a dictionary to hold the columns to add to the general stats table
        headers = OrderedDict()
        # available columns are:
        #SEQ_ID RG  CHIP_ID #SNPS   #READS  AVG_DP  FREEMIX FREELK1 FREELK0 FREE_RH FREE_RA CHIPMIX CHIPLK1 CHIPLK0 CHIP_RH CHIP_RA DPREF   RDPHET  RDPALT
        #,see https://genome.sph.umich.edu/wiki/VerifyBamID#Interpreting_output_files

        # add the CHIPMIX column. set the title and description
        headers['CHIPMIX'] = {
            'title': 'Contamination Prediction - sequence and array',
            'description': 'VerifyBamID: CHIPMIX -   Sequence+array estimate of contamination (NA if the external genotype is unavailable) (0-1 scale)',
            'format': '{:,.5f}' 
        }

        # add the FREEMIX column. set the title and description
        headers['FREEMIX'] = {
            'title': 'Contamination Prediction - sequence only',
            'description': 'VerifyBamID: FREEMIX -   Sequence-only estimate of contamination (0-1 scale).',
            'format': '{:,.5f}' 
        }
        
        # pass the data dictionary and header dictionary to function to add to table.
        self.general_stats_addcols(self.verifybamid_data, headers)

        
    def verifybamid_table(self):
        """
        Create a table with all the columns from verify BAM ID
        """
        # create empty dictionary
        data = dict()
        # for each smaple name and the corresponding dictionary
        for s_name, d in self.verifybamid_data.items():
            # create a entry in data dictionary where key is sample namd and value an empty dictionary
            data[s_name]={}
            # for each column for this sample
            for column in d:
                # add the column name to the sample dictionary and the associated value
                data[s_name][column] = d[column]

        # create an ordered dictionary to preserve the order of columns
        headers = OrderedDict()
        # add each column and the title and description (taken from verifyBAMID website)
        headers['#SEQ_ID']={
            'title': 'Sample Name',
            'description': 'Sample ID of the sequenced sample. Obtained from @RG header / SM tag in the BAM file',
        }
        headers['RG']={
            'title': 'Read Group',
            'description': 'ReadGroup ID of sequenced lane. For [outPrefix].selfSM and [outPrefix].bestSM, these values are "ALL',
        }
        headers['CHIP_ID'] = {
            'title': 'Chip ID',
            'description': 'ReadGroup ID of sequenced lane. For [outPrefix].selfSM and [outPrefix].bestSM, these values are "ALL',
        }
        headers['#SNPS']= {
            'title': '#SNPS',
            'description': '# of SNPs passing the criteria from the VCF file',
        }
        headers['#READS']= {
            'title': '#Reads',
            'description': 'Total # of reads loaded from the BAM file',
        }
        headers['AVG_DP']= {
            'title': 'Average Depth',
            'description': 'Average sequencing depth at the sites in the VCF file',
        }
        # specify the number of decimal points to display
        headers['FREEMIX']= {
            'title': 'FREEMIX',
            'description': 'Sequence-only estimate of contamination (0-1 scale)',
            'format': '{:,.5f}' 
        }
        headers['FREELK1']= {
            'title': 'FREEELK1',
            'description': 'Maximum log-likelihood of the sequence reads given estimated contamination under sequence-only method',
        }
        headers['FREELK0']= {
            'title': 'FREELK0',
            'description': 'Log-likelihood of the sequence reads given no contamination under sequence-only method',
        }
        headers['FREE_RH']= {
            'title': 'FREE_RH',
            'description': 'Estimated reference bias parameter Pr(refBase|HET) (when --free-refBias or --free-full is used)',
        }
        headers['FREE_RA']= {
            'title': 'FREE_RA',
            'description': 'Estimated reference bias parameter Pr(refBase|HOMALT) (when --free-refBias or --free-full is used)',
        }
        # specify the number of decimal points to display
        headers['CHIPMIX']= {
            'title': 'CHIPMIX',
            'description': 'Sequence+array estimate of contamination (NA if the external genotype is unavailable) (0-1 scale)',
            'format': '{:,.5f}' 
        }
        headers[ 'CHIPLK1']= {
            'title': 'CHIPLK1',
            'description': 'Maximum log-likelihood of the sequence reads given estimated contamination under sequence+array method (NA if the external genotypes are unavailable)',
        }
        headers['CHIPLK0']= {
            'title': 'CHIPLK0',
            'description': ' Log-likelihood of the sequence reads given no contamination under sequence+array method (NA if the external genotypes are unavailable)',
        }
        headers['CHIP_RH']= {
            'title': 'CHIP_RH',
            'description': 'Estimated reference bias parameter Pr(refBase|HET) (when --chip-refBias or --chip-full is used)',
        }
        headers['CHIP_RA']= {
            'title': 'CHIP_RA',
            'description': 'Estimated reference bias parameter Pr(refBase|HOMALT) (when --chip-refBias or --chip-full is used)',
        }
        headers['DPREF']= {
            'title': 'DPREF',
            'description': 'Depth (Coverage) of HomRef site (based on the genotypes of (SELF_SM/BEST_SM), passing mapQ, baseQual, maxDepth thresholds.',
        }
        headers['RDPHET']= {
            'title': 'RDPHET',
            'description': 'DPHET/DPREF, Relative depth to HomRef site at Heterozygous site.',
        }
        headers['RDPALT'] = {
            'title': 'RDPALT',
            'description': 'DPHET/DPREF, Relative depth to HomRef site at HomAlt site.',
        }

        # if data was found
        if len(data) > 0:
            # senf the plot to add section function with data dict and headers
            self.add_section (
                anchor = 'verifybamid-table',
                plot = table.plot(data,headers)
            )


	
