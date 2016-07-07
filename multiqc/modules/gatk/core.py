# -*- coding: utf-8 -*-
from collections import OrderedDict
import csv

from multiqc import config, BaseMultiqcModule, plots


class MultiqcModule(BaseMultiqcModule):

    """Parse out information from GATK Variant Eval."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='GATK',
            anchor='gatk',
            target='GATK',
            href='https://www.broadinstitute.org/gatk/',
            info=(" is a toolkit offering a wide variety of tools with a "
                  "primary focus on variant discovery and genotyping."))

        # Parse output for all the files
        samples = {}
        for log_file in self.find_log_files(config.sp['gatk']):
            lines = log_file['f'].splitlines()
            compoverlap_output = compoverlap_parse(lines)
            compoverlap_data = compoverlap_values(compoverlap_output)
            samples[log_file['s_name']] = compoverlap_data

        self.sections = [{
            'name': 'GATK CompOverlap',
            'anchor': 'compoverlap',
            'content': compoverlap_table(samples)
        }]
