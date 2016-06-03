#!/usr/bin/env python

""" MultiQC report module. Holds the output from each
module. Is available to subsequent modules. Contains
helper functions to generate markup for report. """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import fnmatch
import io
import json
import mimetypes
import os
import yaml

from multiqc import logger, config, plots

# Treat defaultdict and OrderedDict as normal dicts for YAML output
from yaml.representer import Representer, SafeRepresenter
yaml.add_representer(defaultdict, Representer.represent_dict)
yaml.add_representer(OrderedDict, Representer.represent_dict)
try:
    yaml.add_representer(unicode, SafeRepresenter.represent_unicode)
except NameError:
    pass # Python 3

# Set up global variables shared across modules
general_stats_data = list()
general_stats_headers = list()
general_stats_html = ''
data_sources = defaultdict(lambda:defaultdict(lambda:defaultdict()))
num_hc_plots = 0
num_mpl_plots = 0
saved_raw_data = dict()

# Make a list of files to search
files = list()
def get_filelist():
    
    def add_file(fn, root):
        
        # Check that we don't want to ignore this file
        i_matches = [n for n in config.fn_ignore_files if fnmatch.fnmatch(fn, n)]
        if len(i_matches) > 0:
            logger.debug("Ignoring file as matched an ignore pattern: {}".format(fn))
            return None
    
        # Use mimetypes to exclude binary files where possible
        (ftype, encoding) = mimetypes.guess_type(os.path.join(root, fn))
        if encoding is not None:
            logger.debug("Ignoring file as is encoded: {}".format(fn))
            return None
        if ftype is not None and ftype.startswith('image'):
            if config.report_imgskips:
                logger.debug("Ignoring file as has filetype '{}': {}".format(ftype, fn))
            return None
        
        # Limit search to files under 5MB to avoid 30GB FastQ files etc.
        try:
            filesize = os.path.getsize(os.path.join(root,fn))
        except (IOError, OSError, ValueError, UnicodeDecodeError):
            logger.debug("Couldn't read file when checking filesize: {}".format(fn))
        else:
            if filesize > config.log_filesize_limit:
                logger.debug("Ignoring file as too large: {}".format(fn))
                return None
        
        # Looks good! Remember this file
        files.append({
            'root': root,
            'fn': fn
        })
    
    # Go through the analysis directories
    for path in config.analysis_dir:
        if os.path.isdir(path):            
            for root, dirnames, filenames in os.walk(path, followlinks=True, topdown=True):
                bname = os.path.basename(root)
                # Skip if this directory name matches config.fn_ignore_dirs
                d_matches = [n for n in config.fn_ignore_dirs if fnmatch.fnmatch(bname, n.rstrip(os.sep))]
                if len(d_matches) > 0:
                    logger.debug("Ignoring directory as matched fn_ignore_dirs: {}".format(bname))
                    continue
                
                # Skip if this directory path matches config.fn_ignore_paths
                p_matches = [n for n in config.fn_ignore_paths if fnmatch.fnmatch(root, n.rstrip(os.sep))]
                if len(p_matches) > 0:
                    logger.debug("Ignoring directory as matched fn_ignore_paths: {}".format(root))
                    continue
                
                # Search filenames in this directory
                for fn in filenames:
                    add_file(fn, root)
        
        elif os.path.isfile(path):
            add_file(os.path.basename(path), os.path.dirname(path))


def general_stats_build_html():
    """ Build the general stats HTML, be that a beeswarm plot or a table. """
    global general_stats_html
    pconfig = {
        'id': 'general_stats_table',
        'table_title': 'General Statistics',
        'save_file': True,
        'raw_data_fn':'multiqc_general_stats'
    }
    general_stats_html = plots.table.plot(general_stats_data, general_stats_headers, pconfig)

def data_sources_tofile ():
    fn = 'multiqc_sources.{}'.format(config.data_format_extensions[config.data_format])
    with io.open (os.path.join(config.data_dir, fn), 'w', encoding='utf-8') as f:
        if config.data_format == 'json':
            jsonstr = json.dumps(data_sources, indent=4, ensure_ascii=False)
            print( jsonstr.encode('utf-8', 'ignore').decode('utf-8'), file=f)
        elif config.data_format == 'yaml':
            yaml.dump(data_sources, f, default_flow_style=False)
        else:
            lines = [['Module', 'Section', 'Sample Name', 'Source']]
            for mod in data_sources:
                for sec in data_sources[mod]:
                    for s_name, source in data_sources[mod][sec].items():
                        lines.append([mod, sec, s_name, source])
            body = '\n'.join(["\t".join(l) for l in lines])
            print( body.encode('utf-8', 'ignore').decode('utf-8'), file=f)
    
    