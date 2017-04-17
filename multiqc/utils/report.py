#!/usr/bin/env python

""" MultiQC report module. Holds the output from each
module. Is available to subsequent modules. Contains
helper functions to generate markup for report. """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import click
import fnmatch
import io
import json
import mimetypes
import os
import re
import yaml

from multiqc import config
logger = config.logger

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

# Make a dict of discovered files for each seach key
files = dict()
def get_filelist():

    # Prep search patterns
    spatterns = [{},{},{},{}]
    for key, sps in config.sp.items():
        files[key] = list()
        if not isinstance(sps, list):
            sps = [sps]
        # Split search patterns according to speed of execution.
        if any([x for x in sps if 'contents' in x]):
            if any([x for x in sps if 'num_lines' in x]):
                spatterns[1][key] = sps
            elif any([x for x in sps if 'max_filesize' in x]):
                spatterns[2][key] = sps
            else:
                spatterns[3][key] = sps
        else:
            spatterns[0][key] = sps

    def add_file(fn, root):
        """
        Function applied to each file found when walking the analysis
        directories. Runs through all search patterns and returns True
        if a match is found.
        """
        f = {'fn': fn, 'root': root}

        # Check that this is a file and not a pipe or anything weird
        if not os.path.isfile(os.path.join(root, fn)):
            return None

        # Check that we don't want to ignore this file
        i_matches = [n for n in config.fn_ignore_files if fnmatch.fnmatch(fn, n)]
        if len(i_matches) > 0:
            logger.debug("Ignoring file as matched an ignore pattern: {}".format(fn))
            return None

        # Limit search to small files, to avoid 30GB FastQ files etc.
        try:
            f['filesize'] = os.path.getsize(os.path.join(root,fn))
        except (IOError, OSError, ValueError, UnicodeDecodeError):
            logger.debug("Couldn't read file when checking filesize: {}".format(fn))
        else:
            if f['filesize'] > config.log_filesize_limit:
                return False

        # Test file for each search pattern
        for patterns in spatterns:
            for key, sps in patterns.items():
                for sp in sps:
                    if search_file (sp, f):
                        # Looks good! Remember this file
                        files[key].append(f)
                        # Don't keep searching this file for other modules
                        if not sp.get('shared', False):
                            return
                        # Don't look at other patterns for this module
                        else:
                            break

    def search_file (pattern, f):
        """
        Function to searach a single file for a single search pattern.
        """
        fn_matched = False
        contents_matched = False

        # Use mimetypes to exclude binary files where possible
        (ftype, encoding) = mimetypes.guess_type(os.path.join(f['root'], f['fn']))
        if encoding is not None:
            return False
        if ftype is not None and ftype.startswith('image'):
            return False

        # Search pattern specific filesize limit
        if pattern.get('max_filesize') is not None and 'filesize' in f:
            if f['filesize'] > pattern.get('max_filesize'):
                return False

        # Search by file name (glob)
        if pattern.get('fn') is not None:
            if fnmatch.fnmatch(f['fn'], pattern['fn']):
                fn_matched = True
                if pattern.get('contents') is None:
                    return True

        # Search by file name (regex)
        if pattern.get('fn_re') is not None:
            if re.match( pattern['fn_re'], f['fn']):
                fn_matched = True
                if pattern.get('contents') is None:
                    return True

        # Search by file contents
        if pattern.get('contents') is not None:
            try:
                with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                    l = 1
                    for line in f:
                        if pattern['contents'] in line:
                            contents_matched = True
                            if pattern.get('fn') is None and pattern.get('fn_re') is None:
                                return True
                            break
                        if pattern.get('num_lines') and l >= pattern.get('num_lines'):
                            break
                        l += 1
            except (IOError, OSError, ValueError, UnicodeDecodeError):
                if config.report_readerrors:
                    logger.debug("Couldn't read file when looking for output: {}".format(fn))
                    return False

        return fn_matched and contents_matched

    # Go through the analysis directories and get file list
    sfiles = []
    for path in config.analysis_dir:
        if os.path.isfile(path):
            sfiles.append([os.path.basename(path), os.path.dirname(path)])
        elif os.path.isdir(path):
            for root, dirnames, filenames in os.walk(path, followlinks=True, topdown=True):
                bname = os.path.basename(root)

                # Skip any sub-directories matching ignore params
                orig_dirnames = dirnames[:]
                for n in config.fn_ignore_dirs:
                    dirnames[:] = [d for d in dirnames if not fnmatch.fnmatch(d, n.rstrip(os.sep))]
                    if len(orig_dirnames) != len(dirnames):
                        removed_dirs = [os.path.join(root, d) for d in set(orig_dirnames).symmetric_difference(set(dirnames))]
                        logger.debug("Ignoring directory as matched fn_ignore_dirs: {}".format(", ".join(removed_dirs)))
                        orig_dirnames = dirnames[:]
                for n in config.fn_ignore_paths:
                    dirnames[:] = [d for d in dirnames if not fnmatch.fnmatch(os.path.join(root, d), n.rstrip(os.sep))]
                    if len(orig_dirnames) != len(dirnames):
                        removed_dirs = [os.path.join(root, d) for d in set(orig_dirnames).symmetric_difference(set(dirnames))]
                        logger.debug("Ignoring directory as matched fn_ignore_paths: {}".format(", ".join(removed_dirs)))

                # Skip *this* directory if matches ignore params
                d_matches = [n for n in config.fn_ignore_dirs if fnmatch.fnmatch(bname, n.rstrip(os.sep))]
                if len(d_matches) > 0:
                    logger.debug("Ignoring directory as matched fn_ignore_dirs: {}".format(bname))
                    continue
                p_matches = [n for n in config.fn_ignore_paths if fnmatch.fnmatch(root, n.rstrip(os.sep))]
                if len(p_matches) > 0:
                    logger.debug("Ignoring directory as matched fn_ignore_paths: {}".format(root))
                    continue
                # Search filenames in this directory
                for fn in filenames:
                    sfiles.append([fn, root])
    # Search through collected files
    with click.progressbar(sfiles, label="Searching {} files..".format(len(sfiles))) as searchfiles:
        for sf in searchfiles:
            add_file(sf[0], sf[1])


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


