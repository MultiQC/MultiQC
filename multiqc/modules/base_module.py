#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

import fnmatch
import io
import logging
import os
import random
import re

from multiqc import plots
from multiqc.utils import report, config, util_functions
logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

class BaseMultiqcModule(object):

    def __init__(self, name='base', anchor='base', target='',href='', info='', extra=''):
        self.name = name
        self.anchor = anchor
        if not target:
            target = self.name
        self.intro = '<p><a href="{0}" target="_blank">{1}</a> {2}</p>{3}'.format(
            href, target, info, extra
        )

    def find_log_files(self, patterns, filecontents=True, filehandles=False):
        """
        Search the analysis directory for log files of interest. Can take either a filename
        suffix or a search string to return only log files that contain relevant info.
        :param patterns: Dict with keys 'fn' or 'contents' (or both). Keys can contain
        string or a list of strings. 'fn' matches filenames, 'contents' matches file contents.
        NB: Both searches return file if *any* of the supplied strings are matched.
        :param filehandles: Set to true to return a file handle instead of slurped file contents
        :return: Yields a set with two items - a sample name generated from the filename
                 and either the file contents or file handle for the current matched file.
                 As yield is used, the results can be iterated over without loading all files at once
        """
        
        # Get the search parameters
        fn_match = None
        contents_match = None
        if 'fn' in patterns:
            fn_match = patterns['fn']
        if 'contents' in patterns:
            contents_match = patterns['contents']
        if fn_match == None and contents_match == None:
            logger.warning("No file patterns specified for find_log_files")
            yield None
                
        # Loop through files, yield results if we find something
        for f in report.files:
            
            # Set up vars
            root = f['root']
            fn = f['fn']
            
            # Make a sample name from the filename
            s_name = self.clean_s_name(fn, root)
            
            # Make search strings into lists if a string is given
            if type(fn_match) is str:
                fn_match = [fn_match]
            if type(contents_match) is str:
                contents_match = [contents_match]
            
            # Search for file names ending in a certain string
            fn_matched = False
            if fn_match is not None:
                for m in fn_match:
                    if fnmatch.fnmatch(fn, m):
                        fn_matched = True
                        if not filehandles and not filecontents:
                            yield {'s_name': s_name, 'root': root, 'fn': fn}
            
            if fn_matched or contents_match is not None:
                try:
                    with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                        
                        # Search this file for our string of interest
                        returnfile = False
                        if contents_match is not None and fn_matched is False:
                            for line in f:
                                for m in contents_match:
                                    if m in line:
                                        returnfile = True
                                        break
                            f.seek(0)
                        else:
                            returnfile = True
                        
                        if returnfile:
                            if filehandles:
                                yield {'s_name': s_name, 'f': f, 'root': root, 'fn': fn}
                            elif filecontents:
                                yield {'s_name': s_name, 'f': f.read(), 'root': root, 'fn': fn}

                except (IOError, OSError, ValueError, UnicodeDecodeError):
                    if config.report_readerrors:
                        logger.debug("Couldn't read file when looking for output: {}".format(fn))
    
    
    def clean_s_name(self, s_name, root):
        """ Helper function to take a long file name and strip it
        back to a clean sample name. Somewhat arbitrary.
        :param s_name: The sample name to clean
        :param root: The directory path that this file is within
        :config.prepend_dirs: boolean, whether to prepend dir name to s_name
        :return: The cleaned sample name, ready to be used
        """
        if root is None:
            root = ''
        # Split then take first section to remove everything after these matches
        for ext in config.fn_clean_exts:
            if type(ext) is str:
                ext = {'type':'truncate', 'pattern':ext}
            if ext['type'] == 'truncate':
                s_name = os.path.basename(s_name.split(ext['pattern'] ,1)[0])
            elif ext['type'] == 'replace':
                s_name = s_name.replace(ext['pattern'], '')
            elif ext['type'] == 'regex':
                re_pattern = re.compile(r'{}'.format(ext['pattern']))
                s_name = re.sub(re_pattern, '', s_name)
            else:
                logger.error('Unrecognised config.fn_clean_exts type: {}'.format(ext['type']))
        if config.prepend_dirs:
            s_name = "{} | {}".format(root.replace(os.sep, ' | '), s_name).lstrip('. | ')
        return s_name
    
    
    def general_stats_addcols(self, data, headers={}, namespace=None):
        """ Helper function to add to the General Statistics variable.
        Adds to report.general_stats and does not return anything. Fills
        in required config variables if not supplied.
        :param data: A dict with the data. First key should be sample name,
                     then the data key, then the data.
        :param headers: Dict / OrderedDict with information for the headers, 
                        such as colour scales, min and max values etc.
                        See docs/writing_python.md for more information.
        :return: None
        """
        # Use the module namespace as the name if not supplied
        if namespace is None:
            namespace = self.name
        
        # Add the module name to the description if not already done
        keys = data.keys()
        if len(headers.keys()) > 0:
            keys = headers.keys()
        for k in keys:
            headers[k]['namespace'] = namespace
            desc = headers[k].get('description', headers[k].get('title', k))
            if 'description' not in headers[k]:
                headers[k]['description'] = desc
        
        # Append to report.general_stats for later assembly into table 
        report.general_stats_data.append(data)
        report.general_stats_headers.append(headers)
    
    def add_data_source(self, f=None, s_name=None, source=None, module=None, section=None):
        try:
            if module is None:
                module = self.name
            if section is None:
                section = 'all_sections'
            if s_name is None:
                s_name = f['s_name']
            if source is None:
                source = os.path.abspath(os.path.join(f['root'], f['fn']))
            report.data_sources[module][section][s_name] = source
        except AttributeError:
            logger.warning('Tried to add data source for {}, but was missing fields data'.format(self.name))
        
    
    def write_data_file(self, data, fn, sort_cols=False, data_format=None):
        """ Saves raw data to a dictionary for downstream use, then redirects
        to report.write_data_file() to create the file in the report directory """
        report.saved_raw_data[fn] = data
        util_functions.write_data_file(data, fn, sort_cols, data_format)
    
    ##################################################
    #### DEPRECIATED FORWARDERS
    def plot_bargraph (self, data, cats=None, pconfig={}):
        """ Depreciated function. Forwards to new location. """
        return plots.bargraph.plot(data, cats, pconfig)
    
    def plot_xy_data(self, data, pconfig={}):
        """ Depreciated function. Forwards to new location. """
        return plots.linegraph.plot(data, pconfig)
