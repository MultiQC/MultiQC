import errno
import logging
import os
import shutil
import sys
import time
import traceback
from importlib.metadata import EntryPoint
from typing import Dict, Union, Callable, List, Optional

import rich
from rich.syntax import Syntax

from multiqc.core.file_search import _required_logs_found
from multiqc.core.utils import _data_tmp_dir, _plots_tmp_dir, _RunError
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.utils import config, report, plugin_hooks, software_versions

logger = logging.getLogger(__name__)


def _run_modules(
    tmp_dir: str,
    run_modules: List[Dict[str, Optional[Dict]]],
    run_module_names: List[str],
    filename: str,
) -> None:
    if filename != "stdout" and config.make_data_dir is True:
        config.data_dir = _data_tmp_dir()
        os.makedirs(config.data_dir)
    else:
        config.data_dir = None
    if filename != "stdout" and config.export_plots is True:
        config.plots_dir = _plots_tmp_dir()
        os.makedirs(config.plots_dir)
    else:
        config.plots_dir = None

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

            # Copy over css & js files if requested by the theme
            try:
                for to, path in report.modules_output[-1].css.items():
                    copy_to = os.path.join(tmp_dir, to)
                    os.makedirs(os.path.dirname(copy_to))
                    shutil.copyfile(path, copy_to)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    pass
                else:
                    raise
            except AttributeError:
                pass
            try:
                for to, path in report.modules_output[-1].js.items():
                    copy_to = os.path.join(tmp_dir, to)
                    os.makedirs(os.path.dirname(copy_to))
                    shutil.copyfile(path, copy_to)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    pass
                else:
                    raise
            except AttributeError:
                pass

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
            shutil.rmtree(tmp_dir)
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
    if not _required_logs_found([m.anchor for m in report.modules_output]):
        raise _RunError()

    # Update report with software versions provided in configs
    software_versions.update_versions_from_config(config, report)

    # Add section for software versions if any are found
    if not config.skip_versions_section and report.software_versions:
        # Importing here to avoid circular imports
        from multiqc.modules.software_versions import MultiqcModule

        report.modules_output.append(MultiqcModule())

    # Special-case module if we want to profile the MultiQC running time
    if config.profile_runtime:
        from multiqc.modules.profile_runtime import MultiqcModule

        report.modules_output.append(MultiqcModule())

    # Did we find anything?
    if len(report.modules_output) == 0:
        logger.warning("No analysis results found. Cleaning up…")
        shutil.rmtree(tmp_dir)
        logger.info("MultiQC complete")
        # Exit with an error code if a module broke
        raise _RunError(sys_exit_code=sys_exit_code)

    if config.make_report:
        # Sort the report module output if we have a config
        if len(getattr(config, "report_section_order", {})) > 0:
            section_id_order = {}
            idx = 10
            for mod in reversed(report.modules_output):
                section_id_order[mod.anchor] = idx
                idx += 10
            for anchor, ss in config.report_section_order.items():
                if anchor not in section_id_order.keys():
                    logger.debug(f"Reordering sections: anchor '{anchor}' not found.")
                    continue
                if ss.get("order") is not None:
                    section_id_order[anchor] = ss["order"]
                if ss.get("after") in section_id_order.keys():
                    section_id_order[anchor] = section_id_order[ss["after"]] + 1
                if ss.get("before") in section_id_order.keys():
                    section_id_order[anchor] = section_id_order[ss["before"]] - 1
            sorted_ids = sorted(section_id_order, key=section_id_order.get)
            report.modules_output = [
                mod for i in reversed(sorted_ids) for mod in report.modules_output if mod.anchor == i
            ]

        # Sort the report sections if we have a config
        # Basically the same as above, but sections within a module
        if len(getattr(config, "report_section_order", {})) > 0:
            # Go through each module
            for midx, mod in enumerate(report.modules_output):
                section_id_order = {}
                # Get a list of the section anchors
                idx = 10
                for s in mod.sections:
                    section_id_order[s["anchor"]] = idx
                    idx += 10
                # Go through each section to be reordered
                for anchor, ss in config.report_section_order.items():
                    # Section to be moved is not in this module
                    if anchor not in section_id_order.keys():
                        logger.debug(f"Reordering sections: anchor '{anchor}' not found for module '{mod.name}'.")
                        continue
                    if ss == "remove":
                        section_id_order[anchor] = False
                        continue
                    if ss.get("order") is not None:
                        section_id_order[anchor] = ss["order"]
                    if ss.get("after") in section_id_order.keys():
                        section_id_order[anchor] = section_id_order[ss["after"]] + 1
                    if ss.get("before") in section_id_order.keys():
                        section_id_order[anchor] = section_id_order[ss["before"]] - 1
                # Remove module sections
                section_id_order = {s: o for s, o in section_id_order.items() if o is not False}
                # Sort the module sections
                sorted_ids = sorted(section_id_order, key=section_id_order.get)
                report.modules_output[midx].sections = [s for i in sorted_ids for s in mod.sections if s["anchor"] == i]

    plugin_hooks.mqc_trigger("after_modules")
