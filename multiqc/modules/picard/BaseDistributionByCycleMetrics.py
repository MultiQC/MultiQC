#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard BaseDistributionByCycleMetrics """

from collections import OrderedDict
import logging
import os
import re

from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)

def read_sample_name(line_iter, clean_fn):
    """
    Consumes lines from the provided line_iter and parses those lines
    as a header for the picard base distribution file.  The header
    file is assumed to contain a line with both 'INPUT' and
    'BaseDistributionByCycle'.

    If the header parses correctly, the sample name is returned.  If
    the header does not parse correctly, None is returned. 
    """
    try:
        while True:
            new_line = line_iter.next()
            new_line = new_line.strip()
            if 'BaseDistributionByCycle' in new_line and 'INPUT' in new_line:
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", new_line)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = clean_fn(s_name)
                    return s_name
    except StopIteration:
        return None

def read_base_distrib_data(line_iter):
    """
    Consumes lines from the provided line_iter and parses those lines
    for base distribution data.  Data should be tab separated and
    immediately preceded by a line of headers:

    READ_END  CYCLE  PCT_A  PCT_C  PCT_G  PCT_T  PCT_N

    Returns either None or a dict mapping cycles to tuples
      (read_end pct_a pct_c pct_g pct_t pct_n)
    where all values are numbers.  

    A None indicates that no lines matching the expected format
    were found.
    """
    try:
        line = line_iter.next()
        while 'BaseDistributionByCycle' not in line and '## METRICS CLASS' not in line:
            line = line_iter.next()
        line = line_iter.next()
        headers = line.strip().split("\t")
        assert headers == ['READ_END', 'CYCLE', 'PCT_A', 'PCT_C', 'PCT_G', 'PCT_T', 'PCT_N']

        # read base distribution by cycle
        data = {}

        row = line_iter.next().strip()
        max_cycle_r1 = None
        while row:
            row_data = list(map(float, row.strip().split("\t")))
            read_end, cycle, pct_a, pct_c, pct_g, pct_t, pct_n = row_data
            cycle = int(cycle)
            if read_end == 1.0:
                if cycle > max_cycle_r1:
                    max_cycle_r1 = cycle
            else:
                cycle = cycle - max_cycle_r1
            data_by_cycle = data.setdefault(read_end, dict())
            data_by_cycle[cycle] = ( 
                pct_a, pct_c, pct_g, pct_t, pct_n
            )
            row = line_iter.next().strip()
        return data
    except StopIteration:
        return None

def parse_reports(self):
    """ Find Picard BaseDistributionByCycleMetrics reports and parse their data """
    
    # Set up vars
    self.picard_baseDistributionByCycle_data = dict()
    self.picard_baseDistributionByCycle_samplestats = dict()
    
    # Go through logs and find Metrics
    base_dist_files = self.find_log_files(
        config.sp['picard']['basedistributionbycycle'],
        filehandles=True
    )

    for f in base_dist_files:
        lines = f['f']

        # read through the header of the file to obtain the
        # sample name
        clean_fn = lambda n: self.clean_s_name(n, f['root'])
        s_name = read_sample_name(lines, clean_fn)
        # pull out the data
        data = read_base_distrib_data(lines)

        # data should be a hierarchical dict
        # data[read_end][cycle]
        assert not (set(data) - set([1, 2]))

        # set up the set of s_names
        if 2 in set(data):
            s_names = {
                1:"%s_R1" % s_name,
                2:"%s_R2" % s_name 
            }
        else:
            s_names = { 1:s_name }

        previously_used = (
            set(s_names.values())&set(self.picard_baseDistributionByCycle_data)
        )

        if previously_used:
            for duped_name in previously_used:
                log.debug(
                    "Duplicate sample name found in {}! "
                    "Overwriting: {}".format(f['fn'], duped_name)
                )
        for name in s_names.values():
            self.add_data_source(f, name, section='BaseDistributionByCycle')

        for read_end in s_names: 
            data_by_cycle = data[read_end]
            s_name = s_names[read_end]
            self.picard_baseDistributionByCycle_data[s_name] = data_by_cycle
            samplestats = {
                'sum_pct_a':0,
                'sum_pct_c':0,
                'sum_pct_g':0,
                'sum_pct_t':0,
                'sum_pct_n':0,
                'cycle_count':0,
            }
            self.picard_baseDistributionByCycle_samplestats[s_name] = samplestats
            for c, row in data_by_cycle.items():
                pct_a, pct_c, pct_g, pct_t, pct_n = row
                samplestats['sum_pct_a'] += pct_a
                samplestats['sum_pct_c'] += pct_c
                samplestats['sum_pct_g'] += pct_g
                samplestats['sum_pct_t'] += pct_t
                samplestats['sum_pct_n'] += pct_n
            samplestats['cycle_count'] += len(data_by_cycle.keys())

    # Calculate summed mean values for all read orientations
    for s_name, v in self.picard_baseDistributionByCycle_samplestats.items():
        v['mean_pct_a'] = v['sum_pct_a'] / v['cycle_count']
        v['mean_pct_c'] = v['sum_pct_c'] / v['cycle_count']
        v['mean_pct_g'] = v['sum_pct_g'] / v['cycle_count']
        v['mean_pct_t'] = v['sum_pct_t'] / v['cycle_count']

    if len(self.picard_baseDistributionByCycle_data) > 0:
        
#        # Write parsed data to a file
#        self.write_data_file(self.picard_baseDistributionByCycle_data, 'multiqc_picard_baseDistributionByCycle')
#        
#        # Do we have median insert sizes?
#        missing_medians = False
#        for v in self.picard_baseDistributionByCycle_samplestats.values():
#            if 'summed_median' not in v:
#                missing_medians = True
        
        # Add to general stats table
        self.general_stats_headers['mean_pct_a'] = {
            'title': 'Mean % A',
            'description': 'Mean % adenine, across all cycles',
            'min': 0,
            'suffix': '%',
            'format': '{:.2f}',
        }
        self.general_stats_headers['mean_pct_g'] = {
            'title': 'Mean % G',
            'description': 'Mean % guanine, across all cycles',
            'min': 0,
            'suffix': '%',
            'format': '{:.2f}',
        }
        self.general_stats_headers['mean_pct_t'] = {
            'title': 'Mean % T',
            'description': 'Mean % thymine, across all cycles',
            'min': 0,
            'suffix': '%',
            'format': '{:.2f}',
        }
        self.general_stats_headers['mean_pct_c'] = {
            'title': 'Mean % C',
            'description': 'Mean % cytosine, across all cycles',
            'min': 0,
            'suffix': '%',
            'format': '{:.2f}',
        }
        for s_name in self.picard_baseDistributionByCycle_samplestats:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.picard_baseDistributionByCycle_samplestats[s_name] )
        
        # Plot the data and add section
        pconfig = {
            'id': 'picard_base_distribution_by_cycle',
            'title': 'Base Distribution',
            'ylab': '%',
            'xlab': 'Cycle #',
            'xDecimals': False,
            'tt_label': '<b>cycle {point.x}</b>: {point.y:.2f} %',
            'ymin': 0,
            'data_labels': [
                {'name': '% adenine', 'ylab': '% adenine'},
                {'name': '% cytosine', 'ylab': '% cytosine'},
                {'name': '% guanine', 'ylab': '% guanine'},
                {'name': '% thymine', 'ylab': '% thymine'},
                {'name': '% undetermined', 'ylab': '% undetermined'},
            ]
        }

        # build list of linegraphs
        linegraph_data = [{}, {}, {}, {}, {}]
        for s_name, cycles in self.picard_baseDistributionByCycle_data.items():
            reformat_items = lambda n: {
                cycle : tup[n] for cycle, tup in cycles.items()
            }
            for linegraph, index in zip(linegraph_data, range(5)):
                linegraph[s_name] = reformat_items(index)

        self.sections.append({
            'name': 'Base Distribution',
            'anchor': 'picard-base-distribution-by-cycle',
            'content': '<p>Plot shows the distribution of bases by cycle.</p>' +
                         plots.linegraph.plot(linegraph_data, pconfig)
        })
    
    
    # Return the number of detected samples to the parent module
    return len(self.picard_baseDistributionByCycle_data)
    
