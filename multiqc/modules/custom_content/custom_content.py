#!/usr/bin/env python

""" Core MultiQC module to parse output from custom script output """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import logging
import json
import os
import yaml

from multiqc import config
from multiqc.utils import report
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table, bargraph, linegraph, scatter, heatmap, beeswarm

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
    cust_mods = defaultdict(lambda: defaultdict(lambda: OrderedDict()))

    # Dictionary to hold search patterns - start with those defined in the config
    search_patterns = ['custom_content']

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
            cust_mods[c_id]['config']['id'] = cust_mods[c_id]['config'].get('id', c_id)
            continue

        # Custom Content ID has search patterns in the config
        if c_id in report.files:
            cust_mods[c_id]['config'] = f
            cust_mods[c_id]['config']['id'] = cust_mods[c_id]['config'].get('id', c_id)
            search_patterns.append(c_id)
            continue

        # We should have had something by now
        log.warn("Found section '{}' in config for under custom_data, but no data or search patterns.".format(c_id))

    # Now go through each of the file search patterns
    bm = BaseMultiqcModule()
    for k in search_patterns:
        for f in bm.find_log_files(k):

            f_extension = os.path.splitext(f['fn'])[1]

            # YAML and JSON files are the easiest
            parsed_data = None
            if f_extension == '.yaml' or f_extension == '.yml':
                try:
                    # Parsing as OrderedDict is slightly messier with YAML
                    # http://stackoverflow.com/a/21048064/713980
                    def dict_constructor(loader, node):
                        return OrderedDict(loader.construct_pairs(node))
                    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)
                    parsed_data = yaml.load(f['f'])
                except Exception as e:
                    log.warning("Error parsing YAML file '{}' (probably invalid YAML)".format(f['fn']))
                    log.warning("YAML error: {}".format(e))
                    break
            elif f_extension == '.json':
                try:
                    # Use OrderedDict for objects so that column order is honoured
                    parsed_data = json.loads(f['f'], object_pairs_hook=OrderedDict)
                except Exception as e:
                    log.warning("Error parsing JSON file '{}' (probably invalid JSON)".format(f['fn']))
                    log.warning("JSON error: {}".format(e))
                    break
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
                    # Update the base config with anything parsed from the file
                    b_config = cust_mods.get(c_id, {}).get('config', {})
                    b_config.update( m_config )
                    # Now set the module config to the merged dict
                    m_config = dict(b_config)
                    s_name = m_config.get('sample_name')
                else:
                    c_id = k
                    m_config = cust_mods.get(c_id, {}).get('config', {})

                # Guess sample name if not given
                if s_name is None:
                    s_name = bm.clean_s_name(f['s_name'], f['root'])

                # Guess c_id if no information known
                if k == 'custom_content':
                    c_id = s_name

                # Add information about the file to the config dict
                if 'files' not in m_config:
                    m_config['files'] = dict()
                m_config['files'].update( { s_name : { 'fn': f['fn'], 'root': f['root'] } } )

                # Guess file format if not given
                if m_config.get('file_format') is None:
                    m_config['file_format'] = _guess_file_format( f )

                # Parse data
                try:
                    parsed_data, conf = _parse_txt( f, m_config )
                    if parsed_data is None or len(parsed_data) == 0:
                        log.warning("Not able to parse custom data in {}".format(f['fn']))
                    else:
                        # Did we get a new section id from the file?
                        if conf.get('id') is not None:
                            c_id = conf.get('id')
                        # heatmap - special data type
                        if type(parsed_data) == list:
                            cust_mods[c_id]['data'] = parsed_data
                        else:
                            cust_mods[c_id]['data'].update(parsed_data)
                        cust_mods[c_id]['config'].update(conf)
                except (IndexError, AttributeError, TypeError):
                    log.error("Unexpected parsing error for {}".format(f['fn']), exc_info=True)
                    raise # testing

    # Filter to strip out ignored sample names
    for k in cust_mods:
        cust_mods[k]['data'] = bm.ignore_samples(cust_mods[k]['data'])

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

        # General Stats
        if mod['config'].get('plot_type') == 'generalstats':
            gsheaders = mod['config'].get('pconfig')
            if gsheaders is None:
                headers = set()
                for d in mod['data'].values():
                    headers.update(d.keys())
                headers = list(headers)
                headers.sort()
                gsheaders = OrderedDict()
                for k in headers:
                    gsheaders[k] = dict()

            # Headers is a list of dicts
            if type(gsheaders) == list:
                hs = OrderedDict()
                for h in gsheaders:
                    for k, v in h.items():
                        hs[k] = v
                gsheaders = hs

            # Add namespace if not specified
            for k in gsheaders:
                if 'namespace' not in gsheaders[k]:
                    gsheaders[k]['namespace'] = c_id

            bm.general_stats_addcols(mod['data'], gsheaders)

        # Initialise this new module class and append to list
        else:
            parsed_modules.append( MultiqcModule(k, mod) )
            log.info("{}: Found {} samples ({})".format(k, len(mod['data']), mod['config'].get('plot_type')))

    # Sort sections if we have a config option for order
    mod_order = getattr(config, 'custom_content', {}).get('order', [])
    sorted_modules = [m for m in parsed_modules if m.anchor not in mod_order ]
    sorted_modules.extend([m for k in mod_order for m in parsed_modules if m.anchor == k ])

    return sorted_modules


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

        pconfig = mod['config'].get('pconfig', {})
        if pconfig.get('title') is None:
            pconfig['title'] = modname

        # Table
        if mod['config'].get('plot_type') == 'table':
            pconfig['sortRows'] = pconfig.get('sortRows', False)
            headers = mod['config'].get('headers')
            self.add_section( plot = table.plot(mod['data'], headers, pconfig) )

        # Bar plot
        elif mod['config'].get('plot_type') == 'bargraph':
            self.add_section( plot = bargraph.plot(mod['data'], mod['config'].get('categories'), pconfig) )

        # Line plot
        elif mod['config'].get('plot_type') == 'linegraph':
            self.add_section( plot = linegraph.plot(mod['data'], pconfig) )

        # Scatter plot
        elif mod['config'].get('plot_type') == 'scatter':
            self.add_section( plot = scatter.plot(mod['data'], pconfig) )

        # Heatmap
        elif mod['config'].get('plot_type') == 'heatmap':
            self.add_section( plot = heatmap.plot(mod['data'], mod['config'].get('xcats'), mod['config'].get('ycats'), pconfig) )

        # Beeswarm plot
        elif mod['config'].get('plot_type') == 'beeswarm':
            self.add_section( plot = beeswarm.plot(mod['data'], pconfig) )

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
    """
    Tries to guess file format, first based on file extension (csv / tsv),
    then by looking for common column separators in the first 10 non-commented lines.
    Splits by tab / comma / space and counts resulting number of columns. Finds the most
    common column count, then comparsed how many lines had this number.
    eg. if tab, all 10 lines should have x columns when split by tab.
    Returns: csv | tsv | spaces   (spaces by default if all else fails)
    """
    filename, file_extension = os.path.splitext(f['fn'])
    tabs = []
    commas = []
    spaces = []
    j = 0
    for l in f['f'].splitlines():
        if not l.startswith('#'):
            j += 1
            tabs.append(len(l.split("\t")))
            commas.append(len(l.split(",")))
            spaces.append(len(l.split()))
        if j == 10:
            break
    tab_mode = max(set(tabs), key=tabs.count)
    commas_mode = max(set(commas), key=commas.count)
    spaces_mode = max(set(spaces), key=spaces.count)
    tab_lc = tabs.count(tab_mode) if tab_mode > 1 else 0
    commas_lc = commas.count(commas_mode) if commas_mode > 1 else 0
    spaces_lc = spaces.count(spaces_mode) if spaces_mode > 1 else 0
    if tab_lc == j:
        return 'tsv'
    elif commas_lc == j:
        return 'csv'
    else:
        if tab_lc > commas_lc and tab_lc > spaces_lc:
            return 'tsv'
        elif commas_lc > tab_lc and commas_lc > spaces_lc:
            return 'csv'
        elif spaces_lc > tab_lc and spaces_lc > commas_lc:
            return 'spaces'
        else:
            if tab_mode == commas_lc and tab_mode > spaces_lc:
                if tab_mode > commas_mode:
                    return 'tsv'
                else:
                    return 'csv'
    return 'spaces'

def _parse_txt(f, conf):
    # Split the data into a list of lists by column
    sep = None
    if conf['file_format'] == 'csv':
        sep = ","
    if conf['file_format'] == 'tsv':
        sep = "\t"
    lines = f['f'].splitlines()
    d = []
    ncols = None
    for l in lines:
        if not l.startswith('#'):
            sections = l.split(sep)
            d.append(sections)
            if ncols is None:
                ncols = len(sections)
            elif ncols != len(sections):
                log.warn("Inconsistent number of columns found in {}! Skipping..".format(f['fn']))
                return (None, conf)

    # Convert values to floats if we can
    first_row_str = 0
    for i,l in enumerate(d):
        for j, v in enumerate(l):
            try:
                d[i][j] = float(v)
            except ValueError:
                if (v.startswith('"') and v.endswith('"')) or (v.startswith("'") and v.endswith("'")):
                    v = v[1:-1]
                d[i][j] = v
                # Count strings in first row (header?)
                if i == 0:
                    first_row_str += 1

    all_numeric = all([ type(l) == float for l in d[i][1:] for i in range(1, len(d)) ])

    # Heatmap: Number of headers == number of lines
    if conf.get('plot_type') is None and first_row_str == len(lines) and all_numeric:
        conf['plot_type'] = 'heatmap'
    if conf.get('plot_type') == 'heatmap':
        conf['xcats'] = d[0][1:]
        conf['ycats'] = [s[0] for s in d[1:]]
        data = [s[1:] for s in d[1:]]
        return (data, conf)

    # Header row of strings, or configured as table
    if first_row_str == len(d[0]) or conf.get('plot_type') == 'table':
        data = OrderedDict()
        for s in d[1:]:
            data[s[0]] = OrderedDict()
            for i, v in enumerate(s[1:]):
                cat = str(d[0][i+1])
                data[s[0]][cat] = v
        # Bar graph or table - if numeric data, go for bar graph
        if conf.get('plot_type') is None:
            allfloats = True
            for r in d:
                for v in r[1:]:
                    allfloats = allfloats and type(v) == float
            if allfloats:
                conf['plot_type'] = 'bargraph'
            else:
                conf['plot_type'] = 'table'
        # Set table col_1 header
        if conf.get('plot_type') == 'table' and d[0][0].strip() != '':
            conf['pconfig'] = conf.get('pconfig', {})
            conf['pconfig']['col1_header'] = d[0][0].strip()
        # Return parsed data
        if conf.get('plot_type') == 'bargraph' or conf.get('plot_type') == 'table':
            return (data, conf)
        else:
            data = OrderedDict() # reset

    # Scatter plot: First row is  str : num : num
    if (conf.get('plot_type') is None and len(d[0]) == 3 and
        type(d[0][0]) != float and type(d[0][1]) == float and type(d[0][2]) == float):
        conf['plot_type'] = 'scatter'

    if conf.get('plot_type') == 'scatter':
        data = dict()
        for s in d:
            try:
                data[s[0]] = {
                    'x': float(s[1]),
                    'y': float(s[2])
                }
            except (IndexError, ValueError):
                pass
        return (data, conf)

    # Single sample line / bar graph - first row has two columns
    if len(d[0]) == 2:
        # Line graph - num : num
        if (conf.get('plot_type') is None and type(d[0][0]) == float and type(d[0][1]) == float):
            conf['plot_type'] = 'linegraph'
        # Bar graph - str : num
        if (conf.get('plot_type') is None and type(d[0][0]) != float and type(d[0][1]) == float):
            conf['plot_type'] = 'bargraph'

        # Data structure is the same
        if (conf.get('plot_type') == 'linegraph' or conf.get('plot_type') == 'bargraph'):
            # Set section id based on directory if not known
            if conf.get('id') is None:
                conf['id'] = os.path.basename(f['root'])
            data = OrderedDict()
            for s in d:
                data[s[0]] = s[1]
            return ( { f['s_name']: data }, conf )

    # Multi-sample line graph: No header row, str : lots of num columns
    if conf.get('plot_type') is None and len(d[0]) > 4 and all_numeric:
        conf['plot_type'] = 'linegraph'

    if conf.get('plot_type') == 'linegraph':
        data = dict()
        # Use 1..n range for x values
        for s in d:
            data[s[0]] = dict()
            for i,v in enumerate(s[1:]):
                j = i+1
                data[s[0]][i+1] = v
        return (data, conf)

    # Got to the end and haven't returned. It's a mystery, capn'!
    log.debug("Not able to figure out a plot type for '{}' ".format(f['fn']) +
      "plot type = {}, all numeric = {}, first row str = {}".format( conf.get('plot_type'), all_numeric, first_row_str ))
    return (None, conf)



