"""
MultiQC lint helpers. Simple additional tests to run when
--strict is specified (outside scope of normal functions)
"""

import logging

from multiqc import report

logger = logging.getLogger(__name__)


def lint_error(msg):
    """Add a lint error to the report"""
    logger.error(msg)
    report.lint_errors.append(msg)
