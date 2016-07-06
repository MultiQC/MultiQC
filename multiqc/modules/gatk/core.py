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


def compoverlap_values(novelties):
    """Parse data per sample for compoverlap output."""
    data = {}
    for novelty in novelties:
        if novelty['Novelty'] == 'all':
            data['reference'] = novelty['CompRod']
            data['comp_rate'] = float(novelty['compRate'])
            data['concordant_rate'] = float(novelty['concordantRate'])
            data['eval_variants'] = int(novelty['nEvalVariants'])
            data['novel_sites'] = int(novelty['novelSites'])
        elif novelty['Novelty'] == 'known':
            data['known_sites'] = int(novelty['nEvalVariants'])
    return data


def compoverlap_parse(lines):
    """Parse out data from CompOverlap tool."""
    all_lines = iter(lines)
    # scan until CompOverlap header
    for line in all_lines:
        if line.startswith('#:GATKTable:CompOverlap'):
            break

    # collect lines until empty lines
    relevant_lines = []
    for raw_line in all_lines:
        line = raw_line.strip()
        if not line:
            break
        else:
            relevant_lines.append(line)

    # parse rows into dict
    values = csv.DictReader(relevant_lines, delimiter=' ', skipinitialspace=True)
    return values


def compoverlap_table(data):
    """Bulid a table from the comp overlaps output."""
    headers = OrderedDict()
    headers['comp_rate'] = {
        'title': 'Comp rate',
        'min': 0,
        'max': 1,
        'format': '{:.2f}',
        'scale': 'Blues',
    }
    headers['concordant_rate'] = {
        'title': 'Concordant rate',
        'min': 0,
        'max': 1,
        'format': '{:.2f}',
        'scale': 'Blues',
    }
    headers['eval_variants'] = {
        'title': 'Evaluated variants',
        'min': 0,
        'format': '{:,}',
    }
    headers['known_sites'] = {
        'title': 'Known sites',
        'min': 0,
        'format': '{:,}',
    }
    headers['novel_sites'] = {
        'title': 'Novel sites',
        'min': 0,
        'format': '{:,}',
    }
    table_config = {
        'id': 'gatk-compoverlap',
        'table_title': 'GTAK CompOverlap',
        'save_file': False,
    }

    table_html = plots.table.plot(data, headers)
    return table_html
