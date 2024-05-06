import logging
from multiqc import config, report

logger = logging.getLogger("multiqc")


class RunResult:
    """
    Returned by a MultiQC run for interactive use. Contains the following information:

    * appropriate error code (e.g. 1 if a module broke, 0 on success)
    * error message if a module broke
    * report instance
    * config instance

    """

    def __init__(self, sys_exit_code: int = 0, message: str = ""):
        self.sys_exit_code = sys_exit_code
        self.message = message
        self.report = report
        self.config = config


class RunError(Exception):
    """
    Used internally in `run` to pass errors from sub-steps.
    """

    def __init__(self, message: str = "", sys_exit_code: int = 1):
        self.message = message
        self.sys_exit_code = sys_exit_code
