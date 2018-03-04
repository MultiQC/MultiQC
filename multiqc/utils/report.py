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
import inspect
import lzstring
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
plot_data = dict()
html_ids = list()
lint_errors = list()
num_hc_plots = 0
num_mpl_plots = 0
saved_raw_data = dict()
last_found_file = None

# Make a dict of discovered files for each seach key
searchfiles = list()
files = dict()
def get_filelist(run_module_names):
    """
    Go through all supplied search directories and assembly a master
    list of files to search. Then fire search functions for each file.
    """
    # Prep search patterns
    spatterns = [{},{},{},{},{},{},{}]
    ignored_patterns = []
    for key, sps in config.sp.items():
        mod_name = key.split('/', 1)[0]
        if mod_name.lower() not in [m.lower() for m in run_module_names]:
            ignored_patterns.append(key)
            continue
        files[key] = list()
        if not isinstance(sps, list):
            sps = [sps]

        # Warn if we have any unrecognised search pattern keys
        unrecognised_keys = [y for x in sps for y in x.keys() if y not in ['fn', 'fn_re', 'contents', 'contents_re', 'num_lines', 'shared', 'max_filesize']]
        if len(unrecognised_keys) > 0:
            logger.warn("Unrecognised search pattern keys for '{}': {}".format(key, ', '.join(unrecognised_keys)))

        # Split search patterns according to speed of execution.
        if any([x for x in sps if 'contents_re' in x]):
            if any([x for x in sps if 'num_lines' in x]):
                spatterns[4][key] = sps
            elif any([x for x in sps if 'max_filesize' in x]):
                spatterns[5][key] = sps
            else:
                spatterns[6][key] = sps
        elif any([x for x in sps if 'contents' in x]):
            if any([x for x in sps if 'num_lines' in x]):
                spatterns[1][key] = sps
            elif any([x for x in sps if 'max_filesize' in x]):
                spatterns[2][key] = sps
            else:
                spatterns[3][key] = sps
        else:
            spatterns[0][key] = sps

    if len(ignored_patterns) > 0:
        logger.debug("Ignored search patterns as didn't match running modules: {}".format(', '.join(ignored_patterns)))

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

    # Go through the analysis directories and get file list
    for path in config.analysis_dir:
        if os.path.isfile(path):
            searchfiles.append([os.path.basename(path), os.path.dirname(path)])
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
                    searchfiles.append([fn, root])
    # Search through collected files
    with click.progressbar(searchfiles, label="Searching {} files..".format(len(searchfiles))) as sfiles:
        for sf in sfiles:
            add_file(sf[0], sf[1])

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
            if pattern.get('contents') is None and pattern.get('contents_re') is None:
                return True

    # Search by file name (regex)
    if pattern.get('fn_re') is not None:
        if re.match( pattern['fn_re'], f['fn']):
            fn_matched = True
            if pattern.get('contents') is None and pattern.get('contents_re') is None:
                return True

    # Search by file contents
    if pattern.get('contents') is not None or pattern.get('contents_re') is not None:
        if pattern.get('contents_re') is not None:
            repattern = re.compile(pattern['contents_re'])
        try:
            with io.open (os.path.join(f['root'],f['fn']), "r", encoding='utf-8') as f:
                l = 1
                for line in f:
                    # Search by file contents (string)
                    if pattern.get('contents') is not None:
                        if pattern['contents'] in line:
                            contents_matched = True
                            if pattern.get('fn') is None and pattern.get('fn_re') is None:
                                return True
                            break
                    # Search by file contents (regex)
                    elif pattern.get('contents_re') is not None:
                        if re.search(repattern, line):
                            contents_matched = True
                            if pattern.get('fn') is None and pattern.get('fn_re') is None:
                                return True
                            break
                    # Break if we've searched enough lines for this pattern
                    if pattern.get('num_lines') and l >= pattern.get('num_lines'):
                        break
                    l += 1
        except (IOError, OSError, ValueError, UnicodeDecodeError):
            if config.report_readerrors:
                logger.debug("Couldn't read file when looking for output: {}".format(f['fn']))
                return False

    return fn_matched and contents_matched

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

def save_htmlid(html_id, skiplint=False):
    """ Take a HTML ID, sanitise for HTML, check for duplicates and save.
    Returns sanitised, unique ID """
    global html_ids
    global lint_errors

    # Trailing whitespace
    html_id_clean = html_id.strip()

    # Trailing underscores
    html_id_clean = html_id_clean.strip('_')

    # Must begin with a letter
    if re.match(r'^[a-zA-Z]', html_id_clean) is None:
        html_id_clean = 'mqc_{}'.format(html_id_clean)

    # Replace illegal characters
    html_id_clean = re.sub('[^a-zA-Z0-9_-]+', '_', html_id_clean)

    # Validate if linting
    if config.lint and not skiplint:
        modname = ''
        codeline = ''
        callstack = inspect.stack()
        for n in callstack:
            if 'multiqc/modules/' in n[1] and 'base_module.py' not in n[1]:
                callpath = n[1].split('multiqc/modules/',1)[-1]
                modname = '>{}< '.format(callpath)
                codeline = n[4][0].strip()
                break
    if config.lint and not skiplint and html_id != html_id_clean:
        errmsg = "LINT: {}HTML ID was not clean ('{}' -> '{}') ## {}".format(modname, html_id, html_id_clean, codeline)
        logger.error(errmsg)
        lint_errors.append(errmsg)

    # Check for duplicates
    i = 1
    html_id_base = html_id_clean
    while html_id_clean in html_ids:
        html_id_clean = '{}-{}'.format(html_id_base, i)
        i += 1
        if config.lint and not skiplint:
            errmsg = "LINT: {}HTML ID was a duplicate ({}) ## {}".format(modname, html_id_clean, codeline)
            logger.error(errmsg)
            lint_errors.append(errmsg)

    # Remember and return
    html_ids.append(html_id_clean)
    return html_id_clean


def compress_json(data):
    """ Take a Python data object. Convert to JSON and compress using lzstring """
    json_string = json.dumps(data).encode('utf-8', 'ignore').decode('utf-8')
    # JSON.parse() doesn't handle `NaN`, but it does handle `null`.
    json_string = json_string.replace('NaN', 'null');
    x = lzstring.LZString()
    return x.compressToBase64(json_string)
