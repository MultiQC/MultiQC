import logging
from multiqc.utils import config
from multiqc.multiqc import run

config.logger = logging.getLogger(__name__)

__version__ = config.version
