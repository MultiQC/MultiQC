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
from pathlib import Path
from typing import Optional

import jinja2

from multiqc import config, report
from multiqc.base_module import Section
from multiqc.core import plugin_hooks, tmp_dir
from multiqc.core.exceptions import NoAnalysisFound
from multiqc.core.log_and_rich import iterate_using_progress_bar
from multiqc.core.tmp_dir import rmtree_with_retries
from multiqc.plots import table
from multiqc.plots.plotly.plot import Plot
from multiqc.utils import megaqc, util_functions

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


def write_results() -> None:
    plugin_hooks.mqc_trigger("before_report_generation")

    # Did we find anything?
    if len(report.modules) == 0:
        raise NoAnalysisFound("No analysis data for any module. Check that input files and directories exist")

    output_names: OutputNames = _set_output_names()

    render_and_export_plots(plots_dir_name=output_names.plots_dir_name)

    if not config.skip_generalstats:
        _render_general_stats_table(plots_dir_name=output_names.plots_dir_name)

    paths: OutputPaths = _create_or_override_dirs(output_names)

    if config.make_data_dir and not paths.to_stdout and paths.data_dir:
        _write_data_files(paths.data_dir)
        logger.info(
            "Data        : {}{}".format(
                _maybe_relative_path(paths.data_dir),
                "   (overwritten)" if paths.data_dir_overwritten else "",
            )
        )
    else:
        logger.info("Data        : None")

    if config.make_report:
        # Render report HTML, write to file or stdout
        _write_html_report(paths.to_stdout, paths.report_path)

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

    if config.export_plots and paths.plots_dir:
        # Copy across the static plot images if requested
        _move_exported_plots(paths.plots_dir)
        logger.info(
            "Plots       : {}{}".format(
                _maybe_relative_path(paths.plots_dir),
                "   (overwritten)" if paths.plots_dir_overwritten else "",
            )
        )

    # Zip the data directory if requested
    if config.zip_data_dir and paths.data_dir is not None:
        shutil.make_archive(str(paths.data_dir), format="zip", root_dir=str(paths.data_dir))
        tmp_dir.rmtree_with_retries(paths.data_dir)

    if paths.report_path:
        logger.debug(f"Report HTML written to {paths.report_path}")


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
        or (config.export_plots and paths.plots_dir and paths.plots_dir.exists())
    ):
        if config.force:
            if config.make_report and isinstance(paths.report_path, Path) and paths.report_path.exists():
                paths.report_overwritten = True
                os.remove(paths.report_path)
            if config.make_data_dir and paths.data_dir and paths.data_dir.exists():
                paths.data_dir_overwritten = True
                shutil.rmtree(paths.data_dir)
            if config.export_plots and paths.plots_dir and paths.plots_dir.exists():
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
            if config.data_dir:
                output_names.data_dir_name = paths.data_dir.name
            if config.plots_dir:
                output_names.plots_dir_name = paths.plots_dir.name
            logger.info("Existing reports found, adding suffix to filenames. Use '--force' to overwrite.")

    if config.make_report and isinstance(paths.report_path, Path):
        paths.report_path.parent.mkdir(exist_ok=True)

    return paths


def render_and_export_plots(plots_dir_name: str):
    """
    Render plot HTML, write PNG/SVG and plot data TSV/JSON to plots_tmp_dir() and data_tmp_dir(). Populates report.plot_data
    """

    def update_fn(_, s: Section):
        if s.plot_id:
            _plot = report.plot_by_id[s.plot_id]
            if isinstance(_plot, Plot):
                s.plot = _plot.add_to_report(plots_dir_name=plots_dir_name)
            elif isinstance(_plot, str):
                s.plot = _plot
            else:
                logger.error(f"Unknown plot type for {s.module}/{s.name}")
        else:
            s.plot = ""

    sections = report.get_all_sections()

    # Show progress bar if writing any flat images, i.e. export_plots requested, or at least one plot is flat
    show_progress = config.export_plots
    if not show_progress:
        for s in sections:
            if s.plot_id:
                plot = report.plot_by_id[s.plot_id]
                if isinstance(plot, Plot) and plot.flat:
                    show_progress = True
                    break

    logger.debug("Rendering plots" + (". This may take a while..." if show_progress else ""))
    iterate_using_progress_bar(
        items=sections,
        update_fn=update_fn,
        item_to_str_fn=lambda s: f"{s.module}/{s.name}" if s.name else s.module,
        desc="rendering plots",
        disable_progress=not show_progress,
    )

    report.some_plots_are_deferred = any(
        isinstance(report.plot_by_id[s.plot_id], Plot) and report.plot_by_id[s.plot_id].defer_render
        for s in sections
        if s.plot_id
    )


def _render_general_stats_table(plots_dir_name: str) -> None:
    """
    Construct HTML for the general stats table.
    """

    # Remove empty data sections from the General Stats table
    empty_keys = [i for i, d in enumerate(report.general_stats_data[:]) if len(d) == 0]
    empty_keys.sort(reverse=True)
    for i in empty_keys:
        del report.general_stats_data[i]
        del report.general_stats_headers[i]

    # Add general-stats IDs to table row headers
    for idx, h in enumerate(report.general_stats_headers):
        for k in h.keys():
            unclean_rid = h[k].get("rid", k)
            rid = re.sub(r"\W+", "_", unclean_rid).strip().strip("_")
            h[k]["rid"] = report.save_htmlid(report.clean_htmlid(rid), skiplint=True)

            ns_html = re.sub(r"\W+", "_", h[k]["namespace"]).strip().strip("_").lower()
            report.general_stats_headers[idx][k]["rid"] = report.save_htmlid(
                f"mqc-generalstats-{ns_html}-{h[k]['rid']}"
            )

    all_hidden = True
    for headers in report.general_stats_headers:
        for h in headers.values():
            if not h.get("hidden", False):
                all_hidden = False
                break

    # Generate the General Statistics HTML & write to file
    if len(report.general_stats_data) > 0 and not all_hidden:
        # Clean previous general stats table if running write_report interactively second time
        if "general_stats_table" in report.html_ids:
            report.html_ids.remove("general_stats_table")
            del report.general_stats_html
        pconfig = {
            "id": "general_stats_table",
            "title": "General Statistics",
            "save_file": True,
            "raw_data_fn": "multiqc_general_stats",
        }
        p = table.plot(report.general_stats_data, report.general_stats_headers, pconfig)
        report.general_stats_html = p.add_to_report(plots_dir_name=plots_dir_name) if isinstance(p, Plot) else p
    else:
        config.skip_generalstats = True


def _write_data_files(data_dir: Path) -> None:
    """
    Write auxiliary data files: module data exports, JSON dump, sources, dev plots data, upload MegaQC
    """

    # Exporting plots to files if requested
    logger.debug("Exporting plot data to files")
    for s in report.get_all_sections():
        if s.plot_id and isinstance(report.plot_by_id.get(s.plot_id), Plot):
            report.plot_by_id[s.plot_id].save_data_files()

    # Modules have run, so data directory should be complete by now. Move its contents.
    logger.debug(f"Moving data file from '{report.data_tmp_dir()}' to '{data_dir}'")

    shutil.copytree(
        report.data_tmp_dir(),
        data_dir,
        # Override default shutil.copy2 function to copy files. The default
        # function copies times and mode, which we want to avoid on purpose
        # to get around the problem with mounted CIFS shares (see #625).
        # shutil.copyfile only copies the file without any metadata.
        copy_function=shutil.copyfile,
    )
    rmtree_with_retries(report.data_tmp_dir())

    # Write the report sources to disk
    report.data_sources_tofile(data_dir)

    # Create a file with the module DOIs
    report.dois_tofile(data_dir, report.modules)

    # Data Export / MegaQC integration - save report data to file or send report data to an API endpoint
    if config.data_dump_file or (config.megaqc_url and config.megaqc_upload):
        dump = report.multiqc_dump_json()
        if config.data_dump_file:
            with (data_dir / "multiqc_data.json").open("w") as f:
                util_functions.dump_json(dump, f, indent=4, ensure_ascii=False)
        if config.megaqc_url:
            megaqc.multiqc_api_post(dump)

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
    rmtree_with_retries(tmp_dir.plots_tmp_dir())


def _write_html_report(to_stdout: bool, report_path: Optional[Path]):
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
        shutil.copytree(parent_template.template_dir, tmp_dir.get_tmp_dir(), dirs_exist_ok=True)

    # Copy the template files to the tmp directory (`dirs_exist_ok` makes sure
    # parent template files are overwritten)
    shutil.copytree(template_mod.template_dir, tmp_dir.get_tmp_dir(), dirs_exist_ok=True)

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


def _write_pdf(report_path: Path) -> Optional[Path]:
    pdf_path = report_path.with_suffix(".pdf")
    pandoc_call = [
        "pandoc",
        "--standalone",
        str(report_path),
        "--output",
        str(pdf_path),
        "--pdf-engine=pdflatex",
        "-V",
        "documentclass=article",
        "-V",
        "geometry=margin=1in",
        "-V",
        "title=",
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
