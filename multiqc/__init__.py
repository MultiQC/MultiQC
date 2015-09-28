import os
import logging
import pkgutil
from pkg_resources import get_distribution

logger = logging.getLogger(__name__)

__version__ = get_distribution("multiqc").version

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

multiqc_submods = [name for _, name, _ in pkgutil.iter_modules(
                    [config.modules_dir])
                ]

# Available templates, in case extras have been added
avail_templates = [ d for d in os.listdir(os.path.join(config.MULTIQC_DIR, 'templates'))
                if os.path.isdir(os.path.join(config.MULTIQC_DIR, 'templates', d)) ]
