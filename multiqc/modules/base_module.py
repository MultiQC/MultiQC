#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

from __future__ import print_function
from collections import OrderedDict
import io
import fnmatch
import logging
import os
import re

from multiqc.utils import report, config, util_functions
logger = logging.getLogger(__name__)

class BaseMultiqcModule(object):

    def __init__(self, name='base', anchor='base', target=None, href=None, info=None, extra=None):
        self.name = name
        self.anchor = anchor
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

    def add_section(self, name=None, anchor=None, description='', helptext='', plot='', content='', autoformat=True):
        """ Add a section to the module report output """

        # Default anchor
        if anchor is None:
            if name is not None:
                nid = name.lower().strip().replace(' ','-')
                anchor = '{}-{}'.format(self.anchor, nid)
            else:
                sl = len(self.sections) + 1
                anchor = '{}-section-{}'.format(self.anchor, sl)

        # Format the content
        if autoformat:
            if len(description) > 0:
                description = '<p class="mqc-section-description">{}</p>'.format(description)
            if len(helptext) > 0:
                helptext = '<p class="mqc-section-helptext">{}</p>'.format(helptext)
            if len(plot) > 0:
                plot = '<div class="mqc-section-plot">{}</div>'.format(plot)

        self.sections.append({
            'name': name,
            'anchor': anchor,
            'description': description,
            'helptext': helptext,
            'plot': plot,
            'content': description + helptext + plot + content
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
                    ext = {'type':'truncate', 'pattern':ext}
                if ext['type'] == 'truncate':
                    s_name = os.path.basename(s_name.split(ext['pattern'] ,1)[0])
                elif ext['type'] == 'replace':
                    s_name = s_name.replace(ext['pattern'], '')
                elif ext['type'] == 'regex':
                    s_name = re.sub(ext['pattern'], '', s_name)
                else:
                    logger.error('Unrecognised config.fn_clean_exts type: {}'.format(ext['type']))
            # Trim off characters at the end of names
            for chrs in config.fn_clean_trim:
                if s_name.endswith(chrs):
                    s_name = s_name[:-len(chrs)]
                if s_name.startswith(chrs):
                    s_name = s_name[len(chrs):]
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
