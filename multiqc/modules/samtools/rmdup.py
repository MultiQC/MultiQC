#!/usr/bin/env python

""" MultiQC module to parse output from Samtools rmdup """

import logging
import re
from collections import OrderedDict
from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)


class RmdupReportMixin():

    def parse_samtools_rmdup(self):
        """ Find Samtools rmdup logs and parse their data """

        self.samtools_rmdup = dict()
        for f in self.find_log_files(config.sp['samtools']['rmdup'], filehandles=True):
            # Example below:
            # [bam_rmdupse_core] 26602816 / 103563641 = 0.2569 in library '   '
            dups_regex = "\[bam_\w+_core\] (\d+) / (\d+) = (\d+\.\d+) in library"
            input_file_regex = "rmdup: Opening (\S+) for read."
            s_name = f['s_name']
            for l in f['f']:
                # try to find name from the log if we can
                match = re.search(input_file_regex, l)
                if match:
                    basefn = os.path.basename(match.group(1))
                    fname, ext = os.path.splitext(basefn)
                    # if it's stdin, then use the filename instead
                    if fname != 'stdin':
                        s_name = self.clean_s_name(fname, f['root'])

                match = re.search(dups_regex, l)
                if match:
                    if s_name in self.samtools_rmdup:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.add_data_source(f, s_name)
                    self.samtools_rmdup[s_name] = dict()
                    self.samtools_rmdup[s_name]['n_dups'] = int(match.group(1))
                    self.samtools_rmdup[s_name]['n_tot'] = int(match.group(2))
                    self.samtools_rmdup[s_name]['n_unique'] = int(match.group(2)) - int(match.group(1))
                    self.samtools_rmdup[s_name]['pct_dups'] = float(match.group(3))*100

        if len(self.samtools_rmdup) > 0:
            # Write parsed report data to a file
            self.write_data_file(self.samtools_rmdup, 'multiqc_samtools_rmdup')
            
            # Make a bar plot showing duplicates
            keys = OrderedDict()
            keys['n_unique'] = {'name': 'Non-duplicated reads'}
            keys['n_dups'] = {'name': 'Duplicated reads'}
            pconfig = {
                'id': 'samtools_rmdup_plot',
                'title': 'Samtools rmdup: Duplicate alignments',
                'yDecimals': False
            }
            self.sections.append({
                'name': 'Duplicates removed',
                'anchor': 'samtools-rmdup',
                'content': plots.bargraph.plot(self.samtools_rmdup, keys, pconfig)
            })
            
            # Add a column to the General Stats table
            # General Stats Table
            stats_headers = OrderedDict()
            stats_headers['pct_dups'] = {
                'title': '% Dups',
                'description': 'Percent of duplicate alignments',
                'min': 0,
                'max': 100,
                'suffix': '%',
                'scale': 'OrRd',
                'format': '{:.1f}%'
            }
            self.general_stats_addcols(self.samtools_rmdup, stats_headers, 'Samtools rmdup')

        return len(self.samtools_rmdup)
