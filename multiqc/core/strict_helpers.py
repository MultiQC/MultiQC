"""
MultiQC lint helpers. Simple additional tests to run when
--strict is specified (outside scope of normal functions)
"""

import logging

from multiqc import config, report

logger = logging.getLogger(__name__)


def lint_error(msg: str):
    if config.strict:
        logger.error(msg)
        report.lint_errors.append(msg)
    else:
        logger.debug(msg)
