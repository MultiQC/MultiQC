import dataclasses

from multiqc.utils import config, report


@dataclasses.dataclass
class RunResult:
    """
    Returned by a MultiQC run for interactive use. Contains the following information:

    * appropriate error code (e.g. 1 if a module broke, 0 on success)
    * error message if a module broke
    * report instance
    * config instance

    """

    sys_exit_code: int = 0
    message: str = ""
    report = report
    config = config


class RunError(Exception):
    """
    Used internally in `run` to pass errors from sub-steps.
    """

    def __init__(self, message: str = "", sys_exit_code: int = 1):
        self.message = message
        self.sys_exit_code = sys_exit_code
