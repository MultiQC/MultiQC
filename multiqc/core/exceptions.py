import logging

logger = logging.getLogger("multiqc")


class RunError(Exception):
    """
    Used internally in `run` to pass errors from sub-steps.
    """

    def __init__(self, message: str = "", sys_exit_code: int = 1):
        self.message = message
        self.sys_exit_code = sys_exit_code


class NoAnalysisFound(Exception):
    """
    Raised when no analysis files were found.
    """

    def __init__(self, message: str = "No analysis results found", sys_exit_code: int = 1):
        self.message = message
        self.sys_exit_code = sys_exit_code
