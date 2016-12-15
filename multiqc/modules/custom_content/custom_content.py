#!/usr/bin/env python

""" Core MultiQC module to parse output from custom script output """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import logging
import json
import os
import yaml

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

def custom_module_classes():
    """
    MultiQC Custom Content class. This module does a lot of different
    things depending on the input and is as flexible as possible.
    NB: THIS IS TOTALLY DIFFERENT TO ALL OTHER MODULES
    """

    # Dict to hold parsed data. Each key should contain a custom data type
    # eg. output from a particular script. Note that this script may pick
    # up many different types of data from many different sources.
    # Second level keys should be 'config' and 'data'. Data key should then
    # contain sample names, and finally data.
    cust_mods = defaultdict(lambda: defaultdict(lambda: dict()))

    # Dictionary to hold search patterns - start with those defined in the config
    search_patterns = OrderedDict()
    search_patterns['core_sp'] = config.sp['custom_content']

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
            cust_mods[c_id]['data'].update( f['data'] )
            cust_mods[c_id]['config'].update( { k:v for k, v in f.items() if k is not 'data' } )
            continue

        # File name patterns supplied in config
        if 'sp' in f:
            cust_mods[c_id]['config'] = f
            search_patterns[c_id] = f['sp']
        else:
            log.debug("Search pattern not found for custom module: {}".format(c_id))

    # Now go through each of the file search patterns
    bm = BaseMultiqcModule()
    for k, sp in search_patterns.items():
        for f in bm.find_log_files(sp):

            f_extension = os.path.splitext(f['fn'])[1]

            # YAML and JSON files are the easiest
            parsed_data = None
            if f_extension == '.yaml' or f_extension == '.yml':
                parsed_data = yaml.load(f['f'])
            elif f_extension == '.json':
                parsed_data = json.loads(f['f'])
            if parsed_data is not None:
                c_id = parsed_data.get('id', k)
                if len(parsed_data.get('data', {})) > 0:
                    cust_mods[c_id]['data'].update( parsed_data['data'] )
                    cust_mods[c_id]['config'].update ( { j:k for j,k in parsed_data.items() if j != 'data' } )
                else:
                    log.warning("No data found in {}".format(f['fn']))

            # txt, csv, tsv etc
            else:
                # Look for configuration details in the header
                m_config = _find_file_header( f )
                s_name = None
                if m_config is not None:
                    c_id = m_config.get('id', k)
                    cust_mods[c_id]['config'].update( m_config )
                    s_name = m_config.get('sample_name')
                else:
                    c_id = k

                # Guess sample name if not given
                if s_name is None:
                    s_name = bm.clean_s_name(f['s_name'], f['root'])

                # Guess c_id if no information known
                if k == 'core_sp':
                    c_id = s_name

                # Add information about the file to the config dict
                if 'files' not in cust_mods[c_id]['config']:
                    cust_mods[c_id]['config']['files'] = dict()
                cust_mods[c_id]['config']['files'].update( { s_name : { 'fn': f['fn'], 'root': f['root'] } } )

                # Guess file format if not given
                if cust_mods[c_id]['config'].get('file_format') is None:
                    cust_mods[c_id]['config']['file_format'] = _guess_file_format( f )

                # Parse data
                parsed_data = _parse_txt( f['f'], cust_mods[c_id]['config'] )
                if parsed_data is None:
                    log.warning("Not able to parse custom data in {}".format(f['fn']))
                else:
                    cust_mods[c_id]['data'][s_name] = parsed_data

                    # Guess plot type if not given
                    if cust_mods[c_id]['config'].get('plot_type') is None:
                        cust_mods[c_id]['config']['plot_type'] = _guess_plot_type( parsed_data )

    # Remove any configs that have no data
    remove_cids = [ k for k in cust_mods if len(cust_mods[k]['data']) == 0 ]
    for k in remove_cids:
        del cust_mods[k]

    if len(cust_mods) == 0:
        log.debug("No custom content found")
        raise UserWarning

    # Go through each data type
    parsed_modules = list()
    for k, mod in cust_mods.items():

        # Initialise this new module class and append to list
        parsed_modules.append( MultiqcModule(k, mod) )
        log.info("{}: Found {} reports".format(k, len(mod['data'])), extra={'mname': k})

    return parsed_modules


class MultiqcModule(BaseMultiqcModule):
    """ Module class, used for each custom content type """

    def __init__(self, c_id, mod):

        modname = mod['config'].get('section_name', c_id.replace('_', ' ').title())

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name = modname,
            anchor = mod['config'].get('section_anchor', c_id),
            href = mod['config'].get('section_href'),
            info = mod['config'].get('description')
        )

        # General Stats
        if mod['config'].get('plot_type') == 'generalstats':
            self.general_stats_addcols(mod['data'], mod['config'].get('pconfig'))

        # Table
        elif mod['config'].get('plot_type') == 'table':
            self.intro += plots.table.plot(mod['data'], mod['config'].get('pconfig'))

        # Bar plot
        elif mod['config'].get('plot_type') == 'bargraph':
            self.intro += plots.bargraph.plot(mod['data'], mod['config'].get('categories'), mod['config'].get('pconfig'))

        # Line plot
        elif mod['config'].get('plot_type') == 'linegraph':
            self.intro += plots.linegraph.plot(mod['data'], mod['config'].get('pconfig'))

        # Scatter plot
        elif mod['config'].get('plot_type') == 'scatter':
            self.intro += plots.scatter.plot(mod['data'], mod['config'].get('pconfig'))

        # Heatmap
        elif mod['config'].get('plot_type') == 'heatmap':
            self.intro += plots.heatmap.plot(mod['data'], mod['config'].get('pconfig'))

        # Beeswarm plot
        elif mod['config'].get('plot_type') == 'beeswarm':
            self.intro += plots.beeswarm.plot(mod['data'], mod['config'].get('pconfig'))

        # Not supplied
        elif mod['config'].get('plot_type') == None:
            log.warning("Plot type not found for content ID '{}'".format(c_id))

        # Not recognised
        else:
            log.warning("Error - custom content plot type '{}' not recognised for content ID {}".format(mod['config'].get('plot_type'), c_id))


def _find_file_header(f):
    # Collect commented out header lines
    hlines = []
    for l in f['f'].splitlines():
        if l.startswith('#'):
            hlines.append(l[1:])
    hconfig = None
    try:
        hconfig = yaml.load("\n".join(hlines))
    except yaml.YAMLError:
        log.debug("Could not parse comment file header for MultiQC custom content: {}".format(f['fn']))
    return hconfig

def _guess_file_format(f):
    filename, file_extension = os.path.splitext(f['fn'])
    if file_extension == 'csv':
        return 'csv'
    if file_extension == 'tsv':
        return 'tsv'
    else:
        # TODO: Look at the file contents to try to guess
        return None

def _guess_plot_type(d):
    col1_str = False
    col2_str = False
    col1_num = False
    col2_num = False
    for k,v in d.items():
        try:
            float(k)
            col1_num = True
        except ValueError:
            col1_str = True
        try:
            float(v)
            col2_num = True
        except ValueError:
            col2_str = True

    if col1_num and col2_num:
        return 'linegraph'
    if col1_str and col2_num:
        return 'bargraph'
    #TODO: Quite a lot of other options to put here.

def _parse_txt(f, conf):
    parsed_data = OrderedDict()
    sep = None
    if conf['file_format'] == 'csv':
        sep = ","
    if conf['file_format'] == 'tsv':
        sep = "\t"
    for l in f.splitlines():
        if not l.startswith('#'):
            s = l.split(sep)
            try:
                try:
                    s0 = float(s[0])
                except ValueError:
                    s0 = s[0]
                try:
                    s1 = float(s[1])
                except ValueError:
                    s1 = s[1]
                parsed_data[s0] = s1
            except IndexError:
                pass
    if len(parsed_data) == 0:
        return None
    else:
        return parsed_data



