#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard MarkDuplicates """

from collections import OrderedDict
import logging
import math
import os
import re

from multiqc import config
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self,
    log_key='picard/markdups',
    section_name='Mark Duplicates',
    section_anchor='picard-markduplicates',
    plot_title='Picard: Deduplication Stats',
    plot_id='picard_deduplication',
    data_filename='multiqc_picard_dups'):
    """ Find Picard MarkDuplicates reports and parse their dataself.
    This function is also used by the biobambam2 module, hence the parameters.
    """

    # Set up vars
    self.picard_dupMetrics_data = dict()

    # Get custom config value
    try:
        merge_multiple_libraries = config.picard_config['markdups_merge_multiple_libraries']
    except (AttributeError, KeyError):
        merge_multiple_libraries = True

    # Go through logs and find Metrics
    for f in self.find_log_files(log_key, filehandles=True):
        s_name = f['s_name']
        base_s_name = f['s_name']
        parsed_data = {}
        keys = None
        in_stats_block = False
        recompute_merged_metrics = False
        for l in f['f']:
            #
            # New log starting
            #
            if 'markduplicates' in l.lower() and 'input' in l.lower():
                # Pull sample name from input
                fn_search = re.search(r"INPUT(?:=|\s+)(\[?[^\s]+\]?)", l, flags=re.IGNORECASE)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1).strip('[]'))
                    s_name = self.clean_s_name(s_name, f['root'])
                    base_s_name = s_name
                continue

            #
            # Start of the METRICS table
            #
            if 'UNPAIRED_READ_DUPLICATES' in l:
                in_stats_block = True
                keys = l.rstrip("\n").split("\t")
                continue

            #
            # Currently parsing the METRICS table
            #
            if in_stats_block:

                # Split the values columns
                vals = l.rstrip("\n").split("\t")

                # End of the METRICS table, or multiple libraries and we're not merging them
                if len(vals) < 6 or (not merge_multiple_libraries and len(parsed_data) > 0):

                    # No data
                    if len(keys) == 0 or len(parsed_data) == 0:
                        break

                    # User has requested each library is kept separate
                    # Update the sample name to append the library name
                    if not merge_multiple_libraries and len(parsed_data) > 0:
                        s_name = '{} - {}'.format(s_name, parsed_data['LIBRARY'])

                    # Skip - No reads
                    try:
                        if parsed_data['READ_PAIRS_EXAMINED'] == 0 and parsed_data['UNPAIRED_READS_EXAMINED'] == 0:
                            log.warn("Skipping MarkDuplicates sample '{}' as log contained no reads".format(s_name))
                            continue
                    # Skip - Missing essential fields
                    except KeyError:
                        log.warn("Skipping MarkDuplicates sample '{}' as missing essential fields".format(s_name))
                        continue

                    # Recompute PERCENT_DUPLICATION and ESTIMATED_LIBRARY_SIZE
                    if recompute_merged_metrics:
                        parsed_data['PERCENT_DUPLICATION'] = calculatePercentageDuplication(parsed_data)
                        parsed_data['ESTIMATED_LIBRARY_SIZE'] = estimateLibrarySize(parsed_data)

                    # Save the data
                    if s_name in self.picard_dupMetrics_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.picard_dupMetrics_data[s_name] = parsed_data
                    self.add_data_source(f, s_name, section='DuplicationMetrics')

                    # End of metrics table - stop parsing
                    if len(vals) < 6:
                        break

                    # On to the next library if not merging
                    else:
                        s_name = base_s_name
                        parsed_data = {}

                #
                # Parse the column values
                #
                for i, k in enumerate(keys):

                    # More than one library present and merging stats
                    if k in parsed_data:
                        recompute_merged_metrics = True
                        try:
                            parsed_data[k] += float(vals[i])
                        except (ValueError, TypeError):
                            parsed_data[k] += " / " + vals[i]

                    # First library
                    else:
                        try:
                            parsed_data[k] = float(vals[i])
                        except ValueError:
                            parsed_data[k] = vals[i]

    #
    # Filter to strip out ignored sample names
    #
    self.picard_dupMetrics_data = self.ignore_samples(self.picard_dupMetrics_data)

    if len(self.picard_dupMetrics_data) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_dupMetrics_data, data_filename)

        # Add to general stats table
        self.general_stats_headers['PERCENT_DUPLICATION'] = {
            'title': '% Dups',
            'description': '{} - Percent Duplication'.format(section_name),
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'OrRd',
            'modify': lambda x: self.multiply_hundred(x)
        }
        for s_name in self.picard_dupMetrics_data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.picard_dupMetrics_data[s_name] )


        # Make the bar plot and add to the MarkDuplicates section
        # NOTE: I had a hard time getting these numbers to add up as expected.
        # If you think I've done something wrong, let me know! Please add an
        # issue here: https://github.com/ewels/MultiQC/issues
        for sn in self.picard_dupMetrics_data.keys():
            self.picard_dupMetrics_data[sn]['UNPAIRED_READ_UNIQUE'] = self.picard_dupMetrics_data[sn]['UNPAIRED_READS_EXAMINED'] - self.picard_dupMetrics_data[sn]['UNPAIRED_READ_DUPLICATES']
            self.picard_dupMetrics_data[sn]['READ_PAIR_NOT_OPTICAL_DUPLICATES'] = self.picard_dupMetrics_data[sn]['READ_PAIR_DUPLICATES'] - self.picard_dupMetrics_data[sn]['READ_PAIR_OPTICAL_DUPLICATES']
            self.picard_dupMetrics_data[sn]['READ_PAIR_UNIQUE'] = self.picard_dupMetrics_data[sn]['READ_PAIRS_EXAMINED'] - self.picard_dupMetrics_data[sn]['READ_PAIR_DUPLICATES']

        keys = OrderedDict()
        keys_r = [
            'READ_PAIR_UNIQUE',
            'UNPAIRED_READ_UNIQUE',
            'READ_PAIR_NOT_OPTICAL_DUPLICATES',
            'READ_PAIR_OPTICAL_DUPLICATES',
            'UNPAIRED_READ_DUPLICATES',
            'UNMAPPED_READS'
        ]
        for k in keys_r:
            keys[k] = {'name': k.replace('_',' ').title()}

        # Config for the plot
        pconfig = {
            'id': plot_id,
            'title': plot_title,
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads',
            'cpswitch_c_active': False
        }

        self.add_section (
            name = section_name,
            anchor = section_anchor,
            plot = bargraph.plot(self.picard_dupMetrics_data, keys, pconfig)
        )

    # Return the number of detected samples to the parent module
    return len(self.picard_dupMetrics_data)


def calculatePercentageDuplication(d):
    """
    Picard calculation to compute percentage duplication

    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob/78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam/DuplicationMetrics.java#L106-L110
    """
    try:
        dups = d['UNPAIRED_READ_DUPLICATES'] + d['READ_PAIR_DUPLICATES'] * 2
        examined = d['UNPAIRED_READS_EXAMINED'] + d['READ_PAIRS_EXAMINED'] * 2
        if examined == 0:
            return 0
        return float(dups) / float(examined)
    except KeyError:
        return 0


def estimateLibrarySize(d):
    """
    Picard calculation to estimate library size

    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob/78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam/DuplicationMetrics.java#L153-L164

    Note: Optical duplicates are contained in duplicates and therefore do not enter the calculation here.
    See also the computation of READ_PAIR_NOT_OPTICAL_DUPLICATES.

     * Estimates the size of a library based on the number of paired end molecules observed
     * and the number of unique pairs observed.
     * <p>
     * Based on the Lander-Waterman equation that states:
     * C/X = 1 - exp( -N/X )
     * where
     * X = number of distinct molecules in library
     * N = number of read pairs
     * C = number of distinct fragments observed in read pairs
    """

    try:
        readPairs = d['READ_PAIRS_EXAMINED'] - d['READ_PAIR_OPTICAL_DUPLICATES']
        uniqueReadPairs = d['READ_PAIRS_EXAMINED'] - d['READ_PAIR_DUPLICATES']
    except KeyError:
        return None

    readPairDuplicates = readPairs - uniqueReadPairs

    if readPairs > 0 and readPairDuplicates > 0:

        m = 1.0;
        M = 100.0;

        if uniqueReadPairs >= readPairs or f(m * uniqueReadPairs, uniqueReadPairs, readPairs) < 0:
            logging.warn("Picard recalculation of ESTIMATED_LIBRARY_SIZE skipped - metrics look wrong")
            return None

        # find value of M, large enough to act as other side for bisection method
        while f(M * uniqueReadPairs, uniqueReadPairs, readPairs) > 0:
            M *= 10.0;

        # use bisection method (no more than 40 times) to find solution
        for i in range(40):
            r = (m + M) / 2.0
            u = f(r * uniqueReadPairs, uniqueReadPairs, readPairs)
            if u == 0:
                break
            elif u > 0:
                m = r
            elif u < 0:
                M = r

        return uniqueReadPairs * (m + M) / 2.0
    else:
        return None

def f(x, c, n):
    """
    Picard calculation used when estimating library size

    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob/78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam/DuplicationMetrics.java#L172-L177

    * Method that is used in the computation of estimated library size.
    """
    return c / x - 1 + math.exp(-n / x)
