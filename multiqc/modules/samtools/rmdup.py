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
        """ Go through log file looking for rmdup output.
        If the
        Grab the name from the RG tag of the preceding bwa command """

        dups_regex = "\[bam_\w+_core\] (\d+) / (\d+) = \((\d+.\d+)\) in library"
        # dups_regex = "rmdup: (Removed|Marked) (\d+) of (\d+) \((\d+.\d+)%\) read ids as duplicates"
        input_file_regex = "rmdup: Opening (\S+) for read."
        # rgtag_name_regex = "\\\\tID:(\S*?)\\\\t"
        data = {}
        s_name = None
        fh = f['f']
        for l in fh:
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
                data['n_dups'] = int(match.group(1))
                data['n_tot'] = int(match.group(2))
                fract_dups = float(match.group(3))
                data['pct_dups'] = fract_dups*100

        if s_name is not None:
            s_name = self.clean_s_name(s_name, f['root'])
            if s_name in self.rmdup_data:
                log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
            self.add_data_source(f, s_name)
            self.rmdup_data[s_name] = data
