import logging
import sys
import time
import traceback
import tracemalloc
from importlib.metadata import EntryPoint
from typing import Dict, Union, Callable, List

import rich
from rich.syntax import Syntax

from multiqc import config, report
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.core import plugin_hooks, software_versions
from multiqc.core.exceptions import RunError, NoAnalysisFound

logger = logging.getLogger(__name__)


def trace_memory(stage: str):
    if config.profile_memory:
        mem_current, mem_peak = tracemalloc.get_traced_memory()
        logger.warning(f"Memory {stage}: {mem_current:,d}b, peak: {mem_peak:,d}b")


def exec_modules(mod_dicts_in_order: List[Dict[str, Dict]]) -> None:
    """
    Execute the modules that have been found and loaded.
    """

    # Only run the modules for which any files were found
    non_empty_modules = {key.split("/")[0].lower() for key, files in report.files.items() if len(files) > 0}
    mod_dicts_in_order = [
        m
        for m in mod_dicts_in_order
        if list(m.keys())[0].lower() in non_empty_modules
        # Always run custom content, as it can have data purely from a MultiQC config file (no search files)
        or list(m.keys())[0].lower() == "custom_content"
    ]
    # Special-case software_versions module must be run in the end
    mod_dicts_in_order = [m for m in mod_dicts_in_order if list(m.keys())[0].lower() != "software_versions"]
    mod_names = [list(m.keys())[0] for m in mod_dicts_in_order]
    if not required_logs_found(mod_names):
        raise RunError()

    # Run the modules!
    plugin_hooks.mqc_trigger("before_modules")
    sys_exit_code = 0
    total_mods_starttime = time.time()

    for mod_idx, mod_dict in enumerate(mod_dicts_in_order):
        mod_starttime = time.time()
        if config.profile_memory:
            tracemalloc.start()

        this_module: str = list(mod_dict.keys())[0]
        logger.debug(f"Running module: {this_module}")
        mod_cust_config: Dict = list(mod_dict.values())[0] or {}
        # noinspection PyBroadException
        try:
            entry_point: EntryPoint = config.avail_modules[this_module]
            module_initializer: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
            setattr(module_initializer, "mod_cust_config", mod_cust_config)
            setattr(module_initializer, "mod_id", this_module)

            # *********************************************
            # RUN MODULE. Heavy part. Run module logic to parse logs and prepare plot data.
            these_modules: Union[BaseMultiqcModule, List[BaseMultiqcModule]] = module_initializer()
            # END RUN MODULE
            # *********************************************

            # Single module initializer can create multiple module objects (see custom_content)
            if not isinstance(these_modules, list):
                these_modules = [these_modules]

            # Clean up non-base attribute to save memory.
            trace_memory("before cleaning up attributes")
            for m in these_modules:
                m.clean_child_attributes()
            trace_memory("after cleaning up attributes")

            # Override duplicated outputs
            for prev_mod in report.modules:
                if prev_mod.name in set(m.name for m in these_modules):
                    logger.info(
                        f'Previous "{prev_mod.name}" run will be overridden. It\'s not yet supported to add new samples to a module with multiqc.parse_logs()'
                    )
                    report.modules.remove(prev_mod)
            report.modules.extend(these_modules)

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
            raise
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive debugging of modules.
                raise

            # Flag the error, but carry on
            class CustomTraceback:
                type, value, traceback = sys.exc_info()

                def __rich_console__(self, console: rich.console.Console, options: rich.console.ConsoleOptions):
                    issue_url = f"https://github.com/MultiQC/MultiQC/issues/new?template=bug_report.md&title={this_module}%20module%20-%20{type.__name__}"
                    err_msg = (
                        f"Please copy this log and report it at [bright_blue][link={issue_url}]"
                        f"https://github.com/MultiQC/MultiQC/issues[/link][/] \n"
                        f"[bold underline]Please attach a file that triggers the error.[/] "
                    )
                    if report.last_found_file:
                        err_msg += f"The last file found was: [green]{report.last_found_file}[/]\n"

                    yield err_msg
                    yield Syntax(traceback.format_exc(), "python")

                def __rich_measure__(self, console: rich.console.Console, options: rich.console.ConsoleOptions):
                    tb_width = max([len(line) for line in traceback.format_exc().split("\n")])
                    log_width = 71
                    if report.last_found_file:
                        log_width += len(report.last_found_file)
                    panel_width = max(tb_width, log_width)
                    return rich.console.Measurement(panel_width, panel_width)

            from multiqc.core.log_and_rich import rich_console_print

            rich_console_print(
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

        report.runtimes.mods[mod_names[mod_idx]] = time.time() - mod_starttime
        if config.profile_memory:
            mem_current, mem_peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            report.peak_memory_bytes_per_module[mod_names[mod_idx]] = mem_peak
            report.diff_memory_bytes_per_module[mod_names[mod_idx]] = mem_current
            logger.warning(
                f"{this_module}: memory change: {mem_current:,d}b, peak during module execution: {mem_peak:,d}b"
            )
        if config.profile_runtime:
            logger.warning(f"{this_module}: module run time: {report.runtimes.mods[mod_names[mod_idx]]:.2f}s")

    report.runtimes.total_mods = time.time() - total_mods_starttime

    # Again, if config.require_logs is set, check if for all explicitly requested
    # modules samples were found.
    if not required_logs_found([m.id for m in report.modules]):
        raise RunError()

    # Update report with software versions provided in configs
    software_versions.update_versions_from_config()

    # Did we find anything?
    if len(report.modules) == 0:
        raise NoAnalysisFound("No analysis results found", sys_exit_code=sys_exit_code)

    plugin_hooks.mqc_trigger("after_modules")


def required_logs_found(modules_with_logs):
    if config.require_logs:
        required_modules_with_no_logs = [
            m
            for m in getattr(config, "run_modules", [])
            if m.lower() not in [m.lower() for m in modules_with_logs]
            and m.lower() not in getattr(config, "exclude_modules", [])
        ]
        if required_modules_with_no_logs:
            logger.critical(
                "The following modules were explicitly requested but no log files were found: {}".format(
                    ", ".join(required_modules_with_no_logs)
                )
            )
            return False
    return True
