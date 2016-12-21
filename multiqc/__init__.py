import logging
from multiqc.utils import config

config.logger = logging.getLogger(__name__)

__version__ = config.version
