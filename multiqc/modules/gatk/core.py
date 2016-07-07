# -*- coding: utf-8 -*-
from multiqc import config, BaseMultiqcModule
from .varianteval import compoverlap, countvariants, titv


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
        for log_file in self.find_log_files(config.sp['gatk']['varianteval']):
            lines = log_file['f'].splitlines()
            compoverlap_output = compoverlap.parse(lines)
            data = compoverlap.values(compoverlap_output)

            count_output = countvariants.parse(lines)
            count_data = countvariants.values(count_output)
            data.update(count_data)

            titv_output = titv.parse(lines)
            titv_data, reference = titv.values(titv_output)
            data.update(titv_data)

            samples[log_file['s_name']] = data

        self.sections = []
        self.sections.append({
            'name': 'Count Variants',
            'anchor': 'countvariants',
            'content': countvariants.plot(samples)
        })
        self.sections.append({
            'name': 'Compare Overlap',
            'anchor': 'compoverlap',
            'content': compoverlap.table(samples)
        })

        general_headers = titv.general_headers(reference)
        self.general_stats_addcols(samples, general_headers)

        self.write_data_file(samples, 'multiqc_gatk')
