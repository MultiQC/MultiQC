import logging
from pkg_resources import get_distribution

logger = logging.getLogger(__name__)

__version__ = get_distribution("multiqc").version

from multiqc.utils import config
from multiqc.modules.base_module import BaseMultiqcModule

config.version = __version__