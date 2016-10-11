#!/usr/bin/env python

""" MultiQC module to parse output from Samtools rmdup """

import logging
import re
from collections import OrderedDict, defaultdict
from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)


class RmdupReportMixin():

    def parse_samtools_rmdup(self):
        """ Find Samtools rmdup logs and parse their data """

        self.samtools_rmdup = dict()
        for f in self.find_log_files(config.sp['samtools']['rmdup']):
            parsed_data = parse_single_report(f['f'])
            if len(parsed_data) > 0:
                if f['s_name'] in self.samtools_rmdup:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f, section='rmdup')
                self.samtools_rmdup[f['s_name']] = parsed_data

        if len(self.samtools_rmdup) > 0:

            # Write parsed report data to a file (restructure first)
            self.write_data_file(self.samtools_rmdup, 'multiqc_samtools_rmdup')

        return len(self.samtools_idxstats)

def parse_single_report(f):
    """ Parse a samtools rmdup """

    """ Go through log file looking for rmdup output.
    If the
    Grab the name from the RG tag of the preceding bwa command """

    # Example below:
    # [bam_rmdupse_core] 26602816 / 103563641 = 0.2569 in library '   '
        
    dups_regex = "\[bam_\w+_core\] (\d+) / (\d+) = \((\d+.\d+)\) in library"
    # dups_regex = "rmdup: (Removed|Marked) (\d+) of (\d+) \((\d+.\d+)%\) read ids as duplicates"
    input_file_regex = "rmdup: Opening (\S+) for read."
    # rgtag_name_regex = "\\\\tID:(\S*?)\\\\t"
    parsed_data = {}
    s_name = None
    for l in f.splitlines():
        # try to find name from RG-tag. If bwa mem is used upstream rmdup with pipes, then the bwa mem command
        # including the read group will be written in the log

        # match = re.search(rgtag_name_regex, l)
        # if match:
        #     s_name = match.group(1)

        # try to find name from the input file name, if used
        match = re.search(input_file_regex, l)
        if match:
            basefn = os.path.basename(match.group(1))
            fname, ext = os.path.splitext(basefn)
            # if it's stdin, then try bwa RG-tag instead
            if fname != 'stdin':
                s_name = fname

        match = re.search(dups_regex, l)
        if match:
            parsed_data['n_dups'] = int(match.group(1))
            parsed_data['n_tot'] = int(match.group(2))
            fract_dups = float(match.group(3))
            parsed_data['pct_dups'] = fract_dups*100

    if s_name is not None:
        s_name = self.clean_s_name(s_name, f['root'])
        if s_name in self.rmdup_data:
            log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
        # self.add_data_source(f, s_name)
        # self.rmdup_data[s_name] = parsed_data
    return parsed_data
