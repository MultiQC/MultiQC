import logging

logger = logging.getLogger(__name__)

from . import config
from . base_module import BaseMultiqcModule


# Order tha modules should appear in report. Try to list in order of analysis,
# eg. FastQC is usually the first step, so should be last in this list.
module_order = [
    # Post-alignment analysis results
    'qualimap', 'featureCounts', 'picard',
    # Alignment tool stats
    'bismark', 'star', 'tophat', 'bowtie2', 'bowtie1',
    # Pre-alignment QC
    'cutadapt', 'fastq_screen', 'fastqc'
]
