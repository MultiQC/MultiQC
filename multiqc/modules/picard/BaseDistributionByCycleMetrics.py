#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard BaseDistributionByCycleMetrics """

import logging
import os
import re

from multiqc.plots import linegraph

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
            new_line = next(line_iter)
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
        line = next(line_iter)
        while 'BaseDistributionByCycle' not in line and '## METRICS CLASS' not in line:
            line = next(line_iter)
        line = next(line_iter)
        headers = line.strip().split("\t")
        assert headers == ['READ_END', 'CYCLE', 'PCT_A', 'PCT_C', 'PCT_G', 'PCT_T', 'PCT_N']

        # read base distribution by cycle
        data = {}

        row = next(line_iter).strip()
        max_cycle_r1 = None
        while row:
            row_data = list(map(float, row.strip().split("\t")))
            read_end, cycle, pct_a, pct_c, pct_g, pct_t, pct_n = row_data
            cycle = int(cycle)
            if read_end == 1.0:
                if max_cycle_r1 is None or cycle > max_cycle_r1:
                    max_cycle_r1 = cycle
            elif max_cycle_r1 is not None:
                cycle = cycle - max_cycle_r1
            data_by_cycle = data.setdefault(read_end, dict())
            data_by_cycle[cycle] = (
                pct_a, pct_c, pct_g, pct_t, pct_n
            )
            row = next(line_iter).strip()
        return data
    except (StopIteration, AssertionError):
        return None

def parse_reports(self):
    """ Find Picard BaseDistributionByCycleMetrics reports and parse their data """

    # Set up vars
    self.picard_baseDistributionByCycle_data = dict()
    self.picard_baseDistributionByCycle_samplestats = dict()

    # Go through logs and find Metrics
    base_dist_files = self.find_log_files('picard/basedistributionbycycle', filehandles=True)

    for f in base_dist_files:
        try:
            lines = iter(f['f'])

            # read through the header of the file to obtain the
            # sample name
            clean_fn = lambda n: self.clean_s_name(n, f['root'])
            s_name = read_sample_name(lines, clean_fn)
            assert s_name is not None

            # pull out the data
            data = read_base_distrib_data(lines)
            assert data is not None

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
        except AssertionError:
            pass

    # Calculate summed mean values for all read orientations
    for s_name, v in self.picard_baseDistributionByCycle_samplestats.items():
        v['mean_pct_a'] = v['sum_pct_a'] / v['cycle_count']
        v['mean_pct_c'] = v['sum_pct_c'] / v['cycle_count']
        v['mean_pct_g'] = v['sum_pct_g'] / v['cycle_count']
        v['mean_pct_t'] = v['sum_pct_t'] / v['cycle_count']

    # Filter to strip out ignored sample names
    self.picard_baseDistributionByCycle_data = self.ignore_samples(self.picard_baseDistributionByCycle_data)

    if len(self.picard_baseDistributionByCycle_data) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_baseDistributionByCycle_samplestats, 'multiqc_picard_baseContent')

        # Plot the data and add section
        pconfig = {
            'id': 'picard_base_distribution_by_cycle',
            'title': 'Base Distribution',
            'ylab': '%',
            'xlab': 'Cycle #',
            'xDecimals': False,
            'tt_label': '<b>cycle {point.x}</b>: {point.y:.2f} %',
            'ymax': 100,
            'ymin': 0,
            'data_labels': [
                {'name': '% Adenine', 'ylab': '% Adenine'},
                {'name': '% Cytosine', 'ylab': '% Cytosine'},
                {'name': '% Guanine', 'ylab': '% Guanine'},
                {'name': '% Thymine', 'ylab': '% Thymine'},
                {'name': '% Undetermined', 'ylab': '% Undetermined'},
            ]
        }

        # build list of linegraphs
        linegraph_data = [{}, {}, {}, {}, {}]
        for s_name, cycles in self.picard_baseDistributionByCycle_data.items():
            reformat_items = lambda n: {
                cycle : tup[n] for cycle, tup in cycles.items()
            }
            for lg, index in zip(linegraph_data, range(5)):
                lg[s_name] = reformat_items(index)

        self.add_section (
            name = 'Base Distribution',
            anchor = 'picard-base-distribution-by-cycle',
            description = 'Plot shows the distribution of bases by cycle.',
            plot = linegraph.plot(linegraph_data, pconfig)
        )


    # Return the number of detected samples to the parent module
    return len(self.picard_baseDistributionByCycle_data)

