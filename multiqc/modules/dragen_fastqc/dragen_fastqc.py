from __future__ import absolute_import

import os

from .base_metrics import DragenBaseMetrics
from .read_metrics import DragenReadMetrics
from .dragen_fastqc_gc_metrics import DragenFastqcGcMetrics
from .content_metrics import DragenContentMetrics
from .util import parse_fastqc_metrics_file

import logging
log = logging.getLogger(__name__)


class MultiqcModule(DragenBaseMetrics, DragenReadMetrics, DragenFastqcGcMetrics, DragenContentMetrics):
    """ DRAGEN provides a number of differrent pipelines and outputs, including base calling, DNA and RNA alignment,
    post-alignment processing and variant calling, covering virtually all stages of typical NGS data processing.
    However, it can be treated as a fast aligner with additional features on top, as users will unlikely use any
    features without enabling DRAGEN mapping. So we will treat this module as an alignment tool module and
    place it accordingly in the module_order list, in docs, etc.

    The QC metrics DRAGEN generates resemble those of samtools-stats, qualimap, mosdepth, bcftools-stats and alike.
    Whenver possible, the visual output is made similar to those modules.

    Note that this MultiQC module supports some of DRAGEN output but not all. Contributions are welcome!

    The code is structured in a way so every mix-in parses one type of QC file that DRAGEN generates
    (e.g. *.mapping_metrics.csv, *.wgs_fine_hist_normal.csv, etc). If a corresponding file is found, a mix-in adds
    a section into the report.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='DRAGEN-FastQc', anchor='DRAGEN-FastQc', target='DRAGEN-FastQc',
            href='https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html',
            info=(" is a Bio-IT Platform that provides ultra-rapid secondary analysis of sequencing data"
                  " using field-programmable gate array technology (FPGA)."))

        self.css = {'assets/css/multiqc_fastqc.css': os.path.join(os.path.dirname(
            __file__), '..', 'fastqc', 'assets', 'css', 'multiqc_fastqc.css')}
        self.js = {
            'assets/js/multiqc_fastqc.js': os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_fastqc.js')}
        self.intro += '<script type="application/json" class="fastqc_passfails">["DRAGEN_FastQc", {"per_base_sequence_content": {"TEST": "pass"}}]</script>'

        data_by_sample = {}
        for f in self.find_log_files('dragen_fastqc'):
            data_by_mate = parse_fastqc_metrics_file(f)
            if f['s_name'] in data_by_sample:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            data_by_sample.update(data_by_mate)

        # Filter to strip out ignored sample names:
        self.fastqc_data = self.ignore_samples(data_by_sample)

        samples_found = set()

        # "POSITIONAL QUALITY" and "POSITIONAL BASE MEAN QUALITY" metrics
        samples_found |= self.add_base_metrics()

        # "READ MEAN QUALITY" and "READ LENGTHS" metrics
        samples_found |= self.add_read_metrics()

        # "GC CONTENT" and "GC CONTENT QUALITY" metrics
        samples_found |= self.add_gc_metrics()

        # "POSITIONAL BASE CONTENT" metrics
        samples_found |= self.add_content_metrics()

        if len(samples_found) == 0:
            raise UserWarning
        log.info("Found {} reports".format(len(samples_found)))
