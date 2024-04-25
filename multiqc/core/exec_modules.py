import logging
import shutil
import sys
import time
import traceback
from importlib.metadata import EntryPoint
from typing import Dict, Union, Callable, List, Optional

import rich
from rich.syntax import Syntax

from multiqc.core.file_search import required_logs_found
from multiqc.core.exceptions import RunError
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.utils import config, report, plugin_hooks, software_versions

logger = logging.getLogger(__name__)


def exec_modules(
    run_modules: List[Dict[str, Optional[Dict]]],
    run_module_names: List[str],
    clean_up: bool = True,
) -> None:
    # Run the modules!
    plugin_hooks.mqc_trigger("before_modules")
    sys_exit_code = 0
    total_mods_starttime = time.time()
    for mod_idx, mod_dict in enumerate(run_modules):
        mod_starttime = time.time()
        this_module: str = list(mod_dict.keys())[0]
        mod_cust_config: Dict = list(mod_dict.values())[0] or {}
        # noinspection PyBroadException
        try:
            entry_point: EntryPoint = config.avail_modules[this_module]
            module_initializer: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
            module_initializer.mod_cust_config = mod_cust_config

            # *********************************************
            # RUN MODULE. Heavy part. Run module logic to parse logs and prepare plot data.
            modules = module_initializer()
            # END RUN MODULE
            # *********************************************

            if not isinstance(modules, list):
                modules = [modules]

            report.modules_output.extend(modules)

        except ModuleNoSamplesFound:
            logger.debug(f"No samples found: {this_module}")
        except UserWarning:  # UserWarning deprecated from 1.16
            msg = f"DEPRECIATED: Please raise 'ModuleNoSamplesFound' instead of 'UserWarning' in module: {this_module}"
            if config.strict:
                logger.error(msg)
                report.lint_errors.append(msg)
            else:
                logger.debug(msg)
            logger.debug(f"No samples found: {this_module}")
        except KeyboardInterrupt:
            if clean_up:
                shutil.rmtree(report.tmp_dir)
            logger.critical(
                "User Cancelled Execution!\n{eq}\n{tb}{eq}\n".format(eq=("=" * 60), tb=traceback.format_exc())
                + "User Cancelled Execution!\nExiting MultiQC..."
            )
            sys.exit(1)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive debugging of modules.
                raise

            # Flag the error, but carry on
            class CustomTraceback:
                def __rich_console__(self, console: rich.console.Console, options: rich.console.ConsoleOptions):
                    sys_tb = sys.exc_info()
                    issue_url = "https://github.com/MultiQC/MultiQC/issues/new?template=bug_report.md&title={}%20module%20-%20{}".format(
                        this_module, sys_tb[0].__name__
                    )
                    yield (
                        "Please copy this log and report it at [bright_blue][link={}]https://github.com/MultiQC/MultiQC/issues[/link][/] \n"
                        "[bold underline]Please attach a file that triggers the error.[/] The last file found was: [green]{}[/]\n".format(
                            issue_url, report.last_found_file
                        )
                    )
                    yield Syntax(traceback.format_exc(), "python")

                def __rich_measure__(self, console: rich.console.Console, options: rich.console.ConsoleOptions):
                    tb_width = max([len(line) for line in traceback.format_exc().split("\n")])
                    try:
                        log_width = 71 + len(report.last_found_file)
                    except TypeError:
                        log_width = 71
                    panel_width = max(tb_width, log_width)
                    return rich.console.Measurement(panel_width, panel_width)

            from multiqc.utils.log import rich_console

            rich_console.print(
                rich.panel.Panel(
                    CustomTraceback(),
                    title=f"Oops! The '[underline]{this_module}[/]' MultiQC module broke...",
                    expand=False,
                    border_style="red",
                    style="on #272822",
                )
            )
            # Still log.debug this so that it ends up in the log file - above is just stderr for now
            logger.debug(
                f"Oops! The '{this_module}' MultiQC module broke...\n"
                + ("=" * 80)
                + "\n"
                + traceback.format_exc()
                + ("=" * 80)
            )
            # Exit code 1 for CI failures etc
            sys_exit_code = 1

        report.runtimes["mods"][run_module_names[mod_idx]] = time.time() - mod_starttime
    report.runtimes["total_mods"] = time.time() - total_mods_starttime

    # Again, if config.require_logs is set, check if for all explicitly requested
    # modules samples were found.
    if not required_logs_found([m.anchor for m in report.modules_output]):
        raise RunError()

    # Update report with software versions provided in configs
    software_versions.update_versions_from_config(config, report)

    # Did we find anything?
    if len(report.modules_output) == 0:
        logger.warning("No analysis results found. Cleaning up…")
        if clean_up:
            shutil.rmtree(report.tmp_dir)
        logger.info("MultiQC complete")
        # Exit with an error code if a module broke
        raise RunError(sys_exit_code=sys_exit_code)

    plugin_hooks.mqc_trigger("after_modules")
