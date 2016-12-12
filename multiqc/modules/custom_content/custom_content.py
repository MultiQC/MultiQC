#!/usr/bin/env python

""" Core MultiQC module to parse output from custom script output """

from __future__ import print_function
import logging
import json
import os
import re
import yaml

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    MultiQC Custom Content class. This module does a lot of different
    things depending on the input and is as flexible as possible.
    """

    def __init__(self):
        
        # Dict to hold parsed data. Each key should contain a custom data type
        # eg. output from a particular script. Note that this script may pick
        # up many different types of data from many different sources.
        # Second level keys should be 'config' and 'data'. Data key should then
        # contain sample names, and finally data.
        self.cust_mods = {}
        
        # Dictionary to hold search patterns - start with those defined in the config
        search_patterns = [ config.sp['custom_content'] ]
        
        # First - find files using patterns described in the config
        config_data = getattr(config, 'custom_data', {})
        for k,f in config_data.items():
            
            # Check that we have a dictionary
            if type(f) != dict:
                log.debug("config.custom_data row was not a dictionary: {}".format(k))
                continue
            c_id = f.get('id', k)
            
            # Data supplied in with config (eg. from a multiqc_config.yaml file in working directory)
            if 'data' in f:
                self.init_data_dict(c_id)
                self.cust_mods[c_id]['data'].update( f['data'] )
                if 'config' in f:
                    self.cust_mods[c_id]['config'].update( f['config'] )
                continue
            
            # File name patterns supplied in config
            sp = {}
            if 'fn' in f:
                sp['fn'] = f['fn']
            if 'contents' in f:
                sp['contents'] = f['contents']
            if len(sp) == 0:
                log.debug("Search pattern not found for custom module: {}".format(c_id))
            else:
                self.init_data_dict(c_id)
                self.cust_mods[c_id]['config'].update(f)
                search_patterns.append(sp)
        
        # Now go through each of the file search patterns
        for sp in search_patterns:
            for f in self.find_log_files(sp):
                
                f_extension = os.path.splitext(f['fn'])[1]
                
                # YAML and JSON files are the easiest
                parsed_data = None
                if f_extension == '.yaml' or f_extension == '.yml':
                    parsed_data = yaml.load(f['f'])
                elif f_extension == '.json':
                    parsed_data = json.loads(f['f'])
                if parsed_data is not None:
                    c_id = parsed_data.get('id', f['fn'])
                    self.init_data_dict(c_id)
                    if len(parsed_data.get('data', {})) > 0:
                        self.cust_mods[c_id]['data'].update( parsed_data['data'] )
                        self.cust_mods[c_id]['config'].update ( { j:k for j,k in parsed_data.items() if j != 'data' } )
                    else:
                        log.warning("No data found in {}".format(f['fn']))
                
                # txt, csv, tsv etc
                else:
                    # Look for configuration details in the header
                    m_config = self.find_file_header( f['f'].splitlines() )
                    log.warning("Parsing custom data that isn't YAML / JSON isn't written yet. Coming soon...")
                    
                    # Guess file format if not given
                    
                    # Guess plot type if not given
                    
                    # Parse data

        
        if len(self.cust_mods) == 0:
            log.debug("No custom content found")
            raise UserWarning
        
        # Go through each data type
        for k, mod in self.cust_mods.items():
            # General Stats
            
            # Table
            
            # Bar plot
            
            # Line plot
            
            # Scatter plot
            
            # Heatmap
            
            # Beeswarm plot
            break
        

        
        print(json.dumps(self.cust_mods, indent=4))

    def init_data_dict(self, c_id):
        if c_id not in self.cust_mods:
            self.cust_mods[c_id] = dict()
        if 'data' not in self.cust_mods[c_id]:
            self.cust_mods[c_id]['data'] = dict()
        if 'config' not in self.cust_mods[c_id]:
            self.cust_mods[c_id]['config'] = dict()
    
    def find_file_header(self, lines):
        for l in lines:
            break
        return None
