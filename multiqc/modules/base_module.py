#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

from __future__ import print_function
from collections import OrderedDict
import io
import fnmatch
import logging
import markdown
import os
import re
import textwrap

from multiqc.utils import report, config, util_functions
logger = logging.getLogger(__name__)

class BaseMultiqcModule(object):

    def __init__(self, name='base', anchor='base', target=None, href=None, info=None, extra=None):

        # Custom options from user config that can overwrite module values
        mod_cust_config = getattr(self, 'mod_cust_config', {})
        self.name = mod_cust_config.get('name', name)
        self.anchor = report.save_htmlid( mod_cust_config.get('anchor', anchor) )
        target = mod_cust_config.get('target', target)
        href = mod_cust_config.get('href', href)
        info = mod_cust_config.get('info', info)
        extra = mod_cust_config.get('extra', extra)

        if info is None:
            info = ''
        if extra is None:
            extra = ''
        if target is None:
            target = self.name
        if href is not None:
            mname = '<a href="{}" target="_blank">{}</a>'.format(href, target)
        else:
            mname = target
        self.intro = '<p>{} {}</p>{}'.format( mname, info, extra )
        self.sections = list()

    def find_log_files(self, sp_key, filecontents=True, filehandles=False):
        """
        Return matches log files of interest.
        :param sp_key: Search pattern key specified in config
        :param filehandles: Set to true to return a file handle instead of slurped file contents
        :return: Yields a dict with filename (fn), root directory (root), cleaned sample name
                 generated from the filename (s_name) and either the file contents or file handle
                 for the current matched file (f).
                 As yield is used, the results can be iterated over without loading all files at once
        """

        # Pick up path filters if specified.
        # Allows modules to be called multiple times with different sets of files
        path_filters = getattr(self, 'mod_cust_config', {}).get('path_filters')

        # Old, depreciated syntax support. Likely to be removed in a future version.
        if isinstance(sp_key, dict):
            report.files[self.name] = list()
            for sf in report.searchfiles:
                if report.search_file(sp_key, {'fn': sf[0], 'root': sf[1]}):
                    report.files[self.name].append({'fn': sf[0], 'root': sf[1]})
            sp_key = self.name
            logwarn = "Depreciation Warning: {} - Please use new style for find_log_files()".format(self.name)
            if len(report.files[self.name]) > 0:
                logger.warn(logwarn)
            else:
                logger.debug(logwarn)
        elif not isinstance(sp_key, str):
            logger.warn("Did not understand find_log_files() search key")
            return

        for f in report.files[sp_key]:

            # If path_filters is given, skip unless match
            if path_filters is not None and len(path_filters) > 0:
                if not all([ fnmatch.fnmatch(f['fn'], pf) for pf in path_filters ]):
                    logger.debug("{} - Skipping '{}' as didn't match module path filters".format(sp_key, f['fn']))
                    continue

            # Make a sample name from the filename
            f['s_name'] = self.clean_s_name(f['fn'], f['root'])
            if filehandles or filecontents:
                try:
                    with io.open (os.path.join(f['root'],f['fn']), "r", encoding='utf-8') as fh:
                        if filehandles:
                            f['f'] = fh
                            yield f
                        elif filecontents:
                            f['f'] = fh.read()
                            yield f
                except (IOError, OSError, ValueError, UnicodeDecodeError):
                    if config.report_readerrors:
                        logger.debug("Couldn't open filehandle when returning file: {}".format(f['fn']))
                        f['f'] = None
            else:
                yield f

    def add_section(self, name=None, anchor=None, description='', helptext='', plot='', content='', autoformat=True, autoformat_type='markdown'):
        """ Add a section to the module report output """

        # Default anchor
        if anchor is None:
            if name is not None:
                nid = name.lower().strip().replace(' ','-')
                anchor = '{}-{}'.format(self.anchor, nid)
            else:
                sl = len(self.sections) + 1
                anchor = '{}-section-{}'.format(self.anchor, sl)

        # Sanitise anchor ID and check for duplicates
        anchor = report.save_htmlid(anchor)

        # Format the content
        if autoformat:
            if len(description) > 0:
                description = textwrap.dedent(description)
                if autoformat_type == 'markdown':
                    description = markdown.markdown(description)
            if len(helptext) > 0:
                helptext = textwrap.dedent(helptext)
                if autoformat_type == 'markdown':
                    helptext = markdown.markdown(helptext)

        # Strip excess whitespace
        description = description.strip()
        helptext = helptext.strip()

        self.sections.append({
            'name': name,
            'anchor': anchor,
            'description': description,
            'helptext': helptext,
            'plot': plot,
            'content': content,
            'print_section': any([ n is not None and len(n) > 0 for n in [description, helptext, plot, content] ])
        })

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
        if config.prepend_dirs:
            sep = config.prepend_dirs_sep
            root = root.lstrip('.{}'.format(os.sep))
            dirs = root.split(os.sep)
            if config.prepend_dirs_depth != 0:
                d_idx = config.prepend_dirs_depth * -1
                if config.prepend_dirs_depth > 0:
                    dirs = dirs[d_idx:]
                else:
                    dirs = dirs[:d_idx]

            s_name = "{}{}{}".format(sep.join(dirs), sep, s_name)
        if config.fn_clean_sample_names:
            # Split then take first section to remove everything after these matches
            for ext in config.fn_clean_exts:
                if type(ext) is str:
                    ext = {'type': 'truncate', 'pattern': ext}
                if ext['type'] == 'truncate':
                    s_name = os.path.basename(s_name.split(ext['pattern'], 1)[0])
                elif ext['type'] in ('remove', 'replace'):
                    if ext['type'] == 'replace':
                        logger.warning("use 'config.fn_clean_sample_names.remove' instead "
                                       "of 'config.fn_clean_sample_names.replace' [deprecated]")
                    s_name = s_name.replace(ext['pattern'], '')
                elif ext['type'] == 'regex':
                    s_name = re.sub(ext['pattern'], '', s_name)
                elif ext['type'] == 'regex_keep':
                    match = re.search(ext['pattern'], s_name)
                    s_name = match.group() if match else s_name
                else:
                    logger.error('Unrecognised config.fn_clean_exts type: {}'.format(ext['type']))
            # Trim off characters at the end of names
            for chrs in config.fn_clean_trim:
                if s_name.endswith(chrs):
                    s_name = s_name[:-len(chrs)]
                if s_name.startswith(chrs):
                    s_name = s_name[len(chrs):]

        # Remove trailing whitespace
        s_name = s_name.strip()

        return s_name

    def ignore_samples(self, data):
        """ Strip out samples which match `sample_names_ignore` """
        try:
            if isinstance(data, OrderedDict):
                newdata = OrderedDict()
            elif isinstance(data, dict):
                newdata = dict()
            else:
                return data
            for k,v in data.items():
                # Match ignore glob patterns
                glob_match = any( fnmatch.fnmatch(k, sn) for sn in config.sample_names_ignore )
                re_match = any( re.match(sn, k) for sn in config.sample_names_ignore_re )
                if not glob_match and not re_match:
                    newdata[k] = v
            return newdata
        except (TypeError, AttributeError):
            return data

    def general_stats_addcols(self, data, headers=None, namespace=None):
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
        if headers is None:
            headers = {}
        # Use the module namespace as the name if not supplied
        if namespace is None:
            namespace = self.name

        # Guess the column headers from the data if not supplied
        if headers is None or len(headers) == 0:
            hs = set()
            for d in data.values():
                hs.update(d.keys())
            hs = list(hs)
            hs.sort()
            headers = OrderedDict()
            for k in hs:
                headers[k] = dict()

        # Add the module name to the description if not already done
        keys = headers.keys()
        for k in keys:
            if 'namespace' not in headers[k]:
                headers[k]['namespace'] = namespace
            if 'description' not in headers[k]:
                headers[k]['description'] = headers[k].get('title', k)

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
    #### DEPRECATED FORWARDERS
    def plot_bargraph (self, data, cats=None, pconfig=None):
        """ Depreciated function. Forwards to new location. """
        from multiqc.plots import bargraph
        if pconfig is None:
            pconfig = {}
        return bargraph.plot(data, cats, pconfig)

    def plot_xy_data(self, data, pconfig=None):
        """ Depreciated function. Forwards to new location. """
        from multiqc.plots import linegraph
        if pconfig is None:
            pconfig = {}
        return linegraph.plot(data, pconfig)
