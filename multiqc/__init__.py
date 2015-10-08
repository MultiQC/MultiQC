import os
import logging
import pkgutil
from pkg_resources import get_distribution

logger = logging.getLogger(__name__)

__version__ = get_distribution("multiqc").version

from . import config
from . base_module import BaseMultiqcModule

config.version = __version__