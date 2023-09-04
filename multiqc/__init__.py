#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" __init__.py
~~~~~~~~~~~~~~~~~~~~
Initialises when multiqc module is loaded.

Makes the following available under the main multiqc namespace:
- run()
- config
- config.logger
- __version__
"""

import logging

from .multiqc import run
from .utils import config

config.logger = logging.getLogger(__name__)

__version__ = config.version
