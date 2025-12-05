import base64
import dataclasses
import errno
import io
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
import uuid
from pathlib import Path
from typing import Optional, cast

import jinja2

from multiqc import config, report
from multiqc.base_module import Section
from multiqc.core import log_and_rich, plot_data_store, plugin_hooks, tmp_dir
from multiqc.core.exceptions import NoAnalysisFound
from multiqc.core.log_and_rich import iterate_using_progress_bar
from multiqc.plots import table
from multiqc.plots.plot import Plot, process_batch_exports
from multiqc.plots.violin import ViolinPlot
from multiqc.types import Anchor
from multiqc.utils import util_functions
from multiqc.utils.util_functions import rmtree_with_retries
from multiqc.utils.material_icons import get_material_icon

logger = logging.getLogger(__name__)


@dataclasses.dataclass
class OutputNames:
    """
    To keep config read-only
    """

    output_fn_name: str
    data_dir_name: str
    plots_dir_name: str


@dataclasses.dataclass
class OutputPaths:
    """
    To keep config read-only
    """

    to_stdout: bool = False

    report_path: Optional[Path] = None
    data_dir: Optional[Path] = None
    plots_dir: Optional[Path] = None

    data_dir_overwritten: bool = False
    plots_dir_overwritten: bool = False
    report_overwritten: bool = False


def write_results(return_html: bool = False) -> Optional[str]:
    plugin_hooks.mqc_trigger("before_report_generation")

    # Did we find anything?
    if len(report.modules) == 0:
        raise NoAnalysisFound("No analysis data for any module. Check that input files and directories exist")

    output_file_names: OutputNames = _set_output_names()

    render_and_export_plots(plots_dir_name=output_file_names.plots_dir_name)

    if not config.skip_generalstats:
        _render_general_stats_table(plots_dir_name=output_file_names.plots_dir_name)

    try:
        report.add_ai_summary()
    except ModuleNotFoundError as e:
        logger.error(e)

    paths: OutputPaths = _create_or_override_dirs(output_file_names)

    if config.make_data_dir and not paths.to_stdout and paths.data_dir:
        _write_data_files(paths.data_dir)
        logger.info(
            "Data        : {}{}".format(
                _maybe_relative_path(paths.data_dir),
                "   (overwritten)" if paths.data_dir_overwritten else "",
            )
        )
    else:
        paths.data_dir = None
        logger.info("Data        : None")

    html_content = None
    if config.make_report:
        # Render report HTML, write to file or stdout
        html_content = _write_html_report(paths.to_stdout, paths.report_path, return_html=return_html)

        if paths.report_path and not config.make_pdf:
            logger.info(
                "Report      : {}{}".format(
                    _maybe_relative_path(paths.report_path),
                    "   (overwritten)" if paths.report_overwritten else "",
                )
            )
        elif paths.report_path and config.make_pdf:
            pdf_path = _write_pdf(paths.report_path)
            if pdf_path:
                logger.info(f"Report      : {_maybe_relative_path(pdf_path)}")

    if paths.plots_dir and tmp_dir.plots_tmp_dir(create=False).exists():
        # Copy across the static plot images if requested
        _move_exported_plots(paths.plots_dir)
        logger.info(
            "Plots       : {}{}".format(
                _maybe_relative_path(paths.plots_dir),
                "   (overwritten)" if paths.plots_dir_overwritten else "",
            )
        )

    # Copy log to the multiqc_data dir. Keeping it in the tmp dir in case if it's an interactive session
    # that goes beyond this write_results run.
    # Do this before zipping the data directory, since zipping will remove the directory.
    if log_and_rich.log_tmp_fn and paths.data_dir and paths.data_dir.exists():
        shutil.copy2(log_and_rich.log_tmp_fn, str(paths.data_dir))

    # Zip the data directory if requested
    if config.zip_data_dir and paths.data_dir is not None:
        shutil.make_archive(str(paths.data_dir), format="zip", root_dir=str(paths.data_dir))
        try:
            util_functions.rmtree_with_retries(paths.data_dir)
        except Exception as e:
            logger.warning(f"Couldn't remove data dir: {e}")

    if paths.report_path:
        logger.debug(f"Report HTML written to {paths.report_path}")

    # Return HTML content if requested
    return html_content if return_html else None


def _maybe_relative_path(path: Path) -> Path:
    """
    If the path is relative to CWD, return the relative path; otherwise, return the full path
    """
    try:
        return path.relative_to(os.getcwd())
    except ValueError:  # could call path.is_relative_to() here, but it's new in Py 3.9
        return path


def _set_output_names() -> OutputNames:
    """
    Set config output paths
    """
    names = OutputNames(
        # These set in config_defaults.yaml or user multiqc_config.yaml, but not through the CLI.
        output_fn_name=config.output_fn_name,
        data_dir_name=config.data_dir_name,
        plots_dir_name=config.plots_dir_name,
    )

    filename = config.filename  # can only be set through the CLI, can be `stdout`
    if filename != "stdout":
        if filename is not None and filename.endswith(".html"):
            filename = filename[:-5]
        if filename is None and config.title is not None:
            filename = re.sub(r"[^\w.-]", "", re.sub(r"[-\s]+", "-", config.title)).strip()
            filename += "_multiqc_report"
        if filename is not None:
            if "output_fn_name" not in config.nondefault_config:
                # Unless the default config.output_fn_name value was explicitly overwritten by user, i.e. set
                # in multiqc_config.yaml distinct to config_defaults.yaml. If not, generating from the command line
                # option `filename`:
                names.output_fn_name = f"{filename}.html"
            if "data_dir_name" not in config.nondefault_config:
                names.data_dir_name = f"{filename}_data"
            if "plots_dir_name" not in config.nondefault_config:
                names.plots_dir_name = f"{filename}_plots"
        if not names.output_fn_name.endswith(".html"):
            names.output_fn_name = f"{names.output_fn_name}.html"

    return names


def _create_or_override_dirs(output_names: OutputNames) -> OutputPaths:
    """
    Write the report data to the output directory
    """

    if config.filename == "stdout":
        logger.info("Printing report to stdout")
        return OutputPaths(to_stdout=True)

    paths = OutputPaths(
        report_path=Path(config.output_fn) if config.output_fn is not None else None,
        data_dir=Path(config.data_dir) if config.data_dir is not None else None,
        plots_dir=Path(config.plots_dir) if config.plots_dir is not None else None,
    )

    output_dir = Path(config.output_dir)

    # Add an output subdirectory if specified by template
    template_mod = config.avail_templates[config.template].load()
    try:
        output_dir = output_dir / template_mod.output_subdir
    except AttributeError:
        pass  # No subdirectory variable given

    if config.make_report:
        paths.report_path = output_dir / output_names.output_fn_name
    else:
        paths.report_path = None
    paths.data_dir = output_dir / output_names.data_dir_name
    paths.plots_dir = output_dir / output_names.plots_dir_name

    # Check for existing reports and remove if -f was specified
    if (
        (config.make_report and isinstance(paths.report_path, Path) and paths.report_path.exists())
        or (config.make_data_dir and paths.data_dir and paths.data_dir.exists())
        or (paths.plots_dir and paths.plots_dir.exists())
    ):
        if config.force:
            if config.make_report and isinstance(paths.report_path, Path) and paths.report_path.exists():
                paths.report_overwritten = True
                os.remove(paths.report_path)
            if config.make_data_dir and paths.data_dir and paths.data_dir.exists():
                paths.data_dir_overwritten = True
                shutil.rmtree(paths.data_dir)
            if paths.plots_dir and paths.plots_dir.exists():
                paths.plots_dir_overwritten = True
                shutil.rmtree(paths.plots_dir)
        else:
            # Set up the base names of the report and the data dir
            report_num = 1

            # Iterate through appended numbers until we find one that's free
            while (
                (config.make_report and isinstance(paths.report_path, Path) and paths.report_path.exists())
                or (config.make_data_dir and paths.data_dir and paths.data_dir.exists())
                or (config.export_plots and paths.plots_dir and paths.plots_dir.exists())
            ):
                if config.make_report:
                    report_base, report_ext = os.path.splitext(output_names.output_fn_name)
                    paths.report_path = output_dir / f"{report_base}_{report_num}{report_ext}"
                if paths.data_dir:
                    dir_base = paths.data_dir.name
                    paths.data_dir = output_dir / f"{dir_base}_{report_num}"
                if paths.plots_dir:
                    plots_base = paths.plots_dir.name
                    paths.plots_dir = output_dir / f"{plots_base}_{report_num}"
                report_num += 1
            if config.make_report and isinstance(paths.report_path, Path):
                output_names.output_fn_name = paths.report_path.name
            if config.data_dir and paths.data_dir:
                output_names.data_dir_name = paths.data_dir.name
            if config.plots_dir and paths.plots_dir:
                output_names.plots_dir_name = paths.plots_dir.name
            logger.info("Existing reports found, adding suffix to filenames. Use '--force' to overwrite.")

    if config.make_report and isinstance(paths.report_path, Path):
        paths.report_path.parent.mkdir(exist_ok=True)

    return paths


def render_and_export_plots(plots_dir_name: str):
    """
    Render plot HTML, write PNG/SVG and plot data TSV/JSON to plots_tmp_dir() and data_tmp_dir().
    Populates report.plot_data
    """

    def update_fn(_, s: Section):
        if s.plot_anchor and s.plot_anchor in report.plot_by_id:
            _plot = report.plot_by_id[s.plot_anchor]
            if isinstance(_plot, Plot):
                s.plot = _plot.add_to_report(
                    plots_dir_name=plots_dir_name,
                    module_anchor=s.module_anchor,
                    section_anchor=s.anchor,
                )
            elif isinstance(_plot, str):
                s.plot = _plot
            else:
                logger.error(f"Unknown plot type for {s.module}/{s.name}")

    sections = report.get_all_sections()

    show_progress = False
    msg = "Rendering plots"
    # Show progress bar if writing any flat images, i.e. export_plots requested, or at least one plot is flat
    if config.export_plots:
        show_progress = True
        msg += (
            f". Export plots to formats {', '.join(config.export_plot_formats)} is requested, so"
            + " it might take a while. To disable plot export, set `export_plots: false` in config,"
            + " or remove the `--export-plots` command line flag"
        )
    elif config.plots_force_flat:
        show_progress = True
        msg += (
            ". Plots are requested to be rendered flat with `--flat` or `plots_force_flat: true` in config,"
            + " and rendering might take a while. To disable, remove or override that flag"
        )
    else:
        for s in sections:
            if s.plot_anchor and s.plot_anchor in report.plot_by_id:
                plot = report.plot_by_id[s.plot_anchor]
                if isinstance(plot, Plot) and plot.flat:
                    show_progress = True
                    msg += (
                        f". Some plots rendered flat because the number of series exceeds {config.plots_flat_numseries},"
                        + " and this rendering might take a while. To disable, set `--interactive`"
                        + " (`plots_force_interactive: true` in config), or change the `plots_flat_numseries`"
                        + " to a higher value"
                    )
                    break

    if show_progress:
        logger.info(msg)
    else:
        logger.debug(msg)

    iterate_using_progress_bar(
        items=sections,
        update_fn=update_fn,
        item_to_str_fn=lambda s: f"{s.module}/{s.name}" if s.name else s.module,
        desc="rendering plots",
        disable_progress=True,
    )

    # Process all batched exports in a single process after all plots are rendered
    process_batch_exports()

    report.some_plots_are_deferred = any(
        isinstance(plot := report.plot_by_id[s.plot_anchor], Plot) and plot.defer_render
        for s in sections
        if s.plot_anchor and s.plot_anchor in report.plot_by_id
    )


def _render_general_stats_table(plots_dir_name: str) -> Optional[Plot]:
    """
    Construct HTML for the general stats table.
    """

    # Remove empty data sections from the General Stats table
    empty_keys = [i for i, d in report.general_stats_data.items() if len(d) == 0]
    empty_keys.sort(reverse=True)
    for i in empty_keys:
        del report.general_stats_data[i]
        del report.general_stats_headers[i]

    # all_hidden = True
    for headers in report.general_stats_headers.values():
        for h in headers.values():
            if not h.get("hidden", False):
                # all_hidden = False
                break

    # Generate the General Statistics HTML & write to file
    # Clean previous general stats table if running write_report interactively second time:
    if Anchor("general_stats_table") in report.html_ids_by_scope[None]:
        report.html_ids_by_scope[None].remove(Anchor("general_stats_table"))  # Violin plot anchor
        if Anchor("general_stats_table_table") in report.html_ids_by_scope[None]:
            report.html_ids_by_scope[None].remove(Anchor("general_stats_table_table"))  # Table anchor
        del report.general_stats_html
    p = table.plot_with_sections(
        data=report.general_stats_data,  # type: ignore
        headers=report.general_stats_headers,  # type: ignore
        pconfig={
            "id": "general_stats_table",
            "title": "General Statistics",
            "save_file": True,
            "raw_data_fn": "multiqc_general_stats",
        },
    )
    if p is None and Anchor("general_stats_table") in report.plot_by_id:
        loaded_plot = report.plot_by_id[Anchor("general_stats_table")]  # loaded from previous run?
        if isinstance(loaded_plot, Plot):
            p = cast(ViolinPlot, loaded_plot)
        elif isinstance(loaded_plot, str):
            p = loaded_plot
        else:
            logger.error("General stats plot is not a Plot object")

    if p is not None:
        if isinstance(p, str):
            report.general_stats_html = p
        else:
            report.plot_by_id[p.anchor] = p
            report.general_stats_html = p.add_to_report(
                plots_dir_name=plots_dir_name,
                module_anchor=Anchor("general_stats_table"),
                section_anchor=Anchor("general_stats_table"),
            )
    else:
        config.skip_generalstats = True
    return None


def _write_data_files(data_dir: Path) -> None:
    """
    Write auxiliary data files: module data exports, JSON dump, sources, dev plots data, upload MegaQC
    """

    # Exporting plots to files if requested
    logger.debug("Exporting plot data to files")
    for s in report.get_all_sections():
        if s.plot_anchor and isinstance(plot := report.plot_by_id.get(s.plot_anchor), Plot):
            plot.save_data_files()

    # Modules have run, so data directory should be complete by now. Move its contents.
    logger.debug(f"Moving data file from '{report.data_tmp_dir()}' to '{data_dir}'")

    # Save metadata to parquet file
    plot_data_store.save_report_metadata()

    shutil.copytree(
        report.data_tmp_dir(),
        data_dir,
        # Override default shutil.copy2 function to copy files. The default
        # function copies times and mode, which we want to avoid on purpose
        # to get around the problem with mounted CIFS shares (see #625).
        # shutil.copyfile only copies the file without any metadata.
        copy_function=shutil.copyfile,
    )
    try:
        rmtree_with_retries(report.data_tmp_dir())
    except Exception as e:
        logger.warning(f"Couldn't remove data tmp dir: {e}")

    # Write the report sources to disk
    report.data_sources_tofile(data_dir)

    # Create a file with the module DOIs
    report.dois_tofile(data_dir, report.modules)

    # Data Export / MegaQC integration - save report data to file or send report data to an API endpoint
    if config.data_dump_file or (config.megaqc_url and config.megaqc_upload):
        report.multiqc_dump_json(data_dir)

    if config.development:
        with (data_dir / "multiqc_plots.js").open("w") as f:
            util_functions.dump_json(report.plot_data, f)


def _move_exported_plots(plots_dir: Path):
    """
    Assuming plots already exported to config.plots_tmp_dir(), move them to config.plots_dir
    """

    # Modules have run, so plots directory should be complete by now. Move its contents.
    logger.debug(f"Moving plots directory from '{tmp_dir.plots_tmp_dir()}' to '{plots_dir}'")

    shutil.copytree(
        tmp_dir.plots_tmp_dir(),
        plots_dir,
        # Override default shutil.copy2 function to copy files. The default
        # function copies times and mode, which we want to avoid on purpose
        # to get around the problem with mounted CIFS shares (see #625).
        # shutil.copyfile only copies the file without any metadata.
        copy_function=shutil.copyfile,
    )
    try:
        rmtree_with_retries(tmp_dir.plots_tmp_dir())
    except Exception as e:
        logger.warning(f"Couldn't remove plots tmp dir: {e}")


def _write_html_report(to_stdout: bool, report_path: Optional[Path], return_html: bool = False) -> Optional[str]:
    """
    Render and write report HTML to disk
    """
    # Copy over css & js files if requested by the theme
    for mod in report.modules:
        if mod.hidden:
            continue
        try:
            for to, path in mod.css.items():
                copy_to = tmp_dir.get_tmp_dir() / to
                copy_to.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, copy_to)
        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise
        except AttributeError:
            pass
        try:
            for to, path in mod.js.items():
                copy_to = tmp_dir.get_tmp_dir() / to
                copy_to.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, copy_to)
        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise
        except AttributeError:
            pass

    plugin_hooks.mqc_trigger("before_template")
    template_mod = config.avail_templates[config.template].load()
    # Load in parent template files first if a child theme
    parent_template = None
    try:
        parent_template = config.avail_templates[template_mod.template_parent].load()
    except AttributeError:
        pass  # Not a child theme
    else:
        shutil.copytree(
            parent_template.template_dir,
            tmp_dir.get_tmp_dir(),
            dirs_exist_ok=True,
            ignore=shutil.ignore_patterns("*.pyc"),
        )

    # Copy the template files to the tmp directory (`dirs_exist_ok` makes sure
    # parent template files are overwritten)
    shutil.copytree(
        template_mod.template_dir, tmp_dir.get_tmp_dir(), dirs_exist_ok=True, ignore=shutil.ignore_patterns("*.pyc")
    )

    # Function to include file contents in Jinja template
    def include_file(name, fdir=tmp_dir.get_tmp_dir(), b64=False):
        try:
            if fdir is None:
                fdir = ""
            _path: str = os.path.join(fdir, name)

            if config.development:
                if os.path.exists(dev_path := os.path.join(template_mod.template_dir, name)):
                    fdir = template_mod.template_dir
                    name = dev_path
                    _path = dev_path
                elif parent_template and os.path.exists(dev_path := os.path.join(parent_template.template_dir, name)):
                    fdir = template_mod.template_dir
                    name = dev_path
                    _path = dev_path

                if re.match(r".*\.min\.(js|css)$", name):
                    unminimized_name = re.sub(r"\.min\.", ".", name)
                    if os.path.exists(os.path.join(fdir, unminimized_name)):
                        name = unminimized_name

                if name.endswith(".js"):
                    return f'</script><script type="text/javascript" src="{name}">'
                if name.endswith(".css"):
                    return f'</style><link rel="stylesheet" href="{name}">'

            if b64:
                with io.open(_path, "rb") as f:
                    return base64.b64encode(f.read()).decode("utf-8")
            else:
                with io.open(_path, "r", encoding="utf-8") as f:
                    return f.read()
        except (OSError, IOError) as e:
            logger.error(f"Could not include file '{name}': {e}")

    # Load the report template
    try:
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(tmp_dir.get_tmp_dir()))
        env.globals["include_file"] = include_file

        # Add Material Design Icons function to all templates
        env.globals["material_icon"] = get_material_icon

        # Add template functions if available
        if hasattr(template_mod, "template_functions"):
            for func_name, func in template_mod.template_functions.items():
                env.globals[func_name] = func

        j_template = env.get_template(template_mod.base_fn, globals={"development": config.development})
    except:  # noqa: E722
        raise IOError(f"Could not load {config.template} template file '{template_mod.base_fn}'")

    # Compress the report plot JSON data
    runtime_compression_start = time.time()
    logger.debug("Compressing plot data")
    report.plot_compressed_json = report.compress_json(report.plot_data)
    report.runtimes.total_compression = time.time() - runtime_compression_start

    # Use jinja2 to render the template and overwrite
    report.analysis_files = [os.path.realpath(d) for d in report.analysis_files]
    report.report_uuid = str(uuid.uuid4())

    # Allow templates to override config settings
    if hasattr(template_mod, "template_dark_mode"):
        config.template_dark_mode = template_mod.template_dark_mode
    if hasattr(template_mod, "plot_font_family"):
        config.plot_font_family = template_mod.plot_font_family

    report_output = j_template.render(report=report, config=config)
    if to_stdout:
        print(report_output, file=sys.stdout)
    else:
        assert report_path is not None
        try:
            with io.open(report_path, "w", encoding="utf-8") as f:
                print(report_output, file=f)
        except IOError as e:
            raise IOError(f"Could not print report to '{config.output_fn}' - {IOError(e)}")

        # Copy over files if requested by the theme
        try:
            for copy_file in template_mod.copy_files:
                fn = tmp_dir.get_tmp_dir() / copy_file
                dest_dir = report_path.parent / copy_file
                shutil.copytree(fn, dest_dir, dirs_exist_ok=True)
        except AttributeError:
            pass  # No files to copy

    # Return HTML content if requested
    return report_output if return_html else None


def _write_pdf(report_path: Path) -> Optional[Path]:
    pdf_path = report_path.with_suffix(".pdf")
    pandoc_call = [
        "pandoc",
        "--standalone",
        str(report_path),
        "--output",
        str(pdf_path),
        "--pdf-engine=lualatex",
        "-V",
        "documentclass=article",
        "-V",
        "geometry=margin=1in",
        "-V",
        "mainfont=DejaVu Sans",
        "-V",
        "sansfont=DejaVu Sans",
        "-V",
        "monofont=DejaVu Sans Mono",
        "-V",
        "fontsize=10pt",
        "-V",
        "title=",
        "-V",
        "tables=true",
    ]
    if config.pandoc_template is not None:
        pandoc_call.append(f"--template={config.pandoc_template}")
    logger.debug(f"Attempting Pandoc conversion to PDF with following command:\n{' '.join(pandoc_call)}")

    try:
        pdf_exit_code = subprocess.call(pandoc_call)
    except OSError as e:
        if e.errno == errno.ENOENT:
            logger.error("Error creating PDF - `pandoc` not found. Is it installed? http://pandoc.org/")
        else:
            logger.error(
                "Error creating PDF! Something went wrong when creating the PDF\n"
                + ("=" * 60)
                + f"\n{traceback.format_exc()}\n"
                + ("=" * 60)
            )
        return None

    if pdf_exit_code != 0:
        logger.error("Error creating PDF! Pandoc returned a non-zero exit code.")
        return None

    # Remove the HTML report
    os.remove(report_path)
    return pdf_path
