import logging

logger = logging.getLogger(__name__)

from multiqc.utils import config
from multiqc.modules.base_module import BaseMultiqcModule

__version__ = config.version
