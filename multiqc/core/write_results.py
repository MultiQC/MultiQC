import base64
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
from typing import Set

from multiqc.core.file_search import include_or_exclude_modules
from multiqc.plots import table

import jinja2

from multiqc import config, report
from multiqc.core.exceptions import RunError
from multiqc.plots.plotly.plot import Plot
from multiqc.core import plugin_hooks
from multiqc.utils import megaqc, util_functions

logger = logging.getLogger(__name__)


def write_results(clean_up=True) -> None:
    plugin_hooks.mqc_trigger("before_report_generation")

    # Did we find anything?
    if len(report.modules) == 0:
        logger.warning("No analysis results found to make a report")
        return

    if config.make_report:
        _order_modules_and_sections()

    _set_output_paths()

    render_and_export_plots()

    if not config.skip_generalstats:
        _render_general_stats_table()

    overwritten = _create_or_override_dirs() if config.filename != "stdout" else set()

    if config.make_data_dir and config.filename != "stdout":
        _write_data_files()
        logger.info(
            "Data        : {}{}".format(
                os.path.relpath(config.data_dir),
                "   (overwritten)" if "data_dir" in overwritten else "",
            )
        )
    else:
        logger.info("Data        : None")

    if config.make_report:
        # Render report HTML, write to file or stdout
        _write_report()
        if config.filename != "stdout":
            if config.make_report:
                logger.info(
                    "Report      : {}{}".format(
                        os.path.relpath(config.output_fn),
                        "   (overwritten)" if "report" in overwritten else "",
                    )
                )
                logger.debug(f"Full report path: {os.path.realpath(config.output_fn)}")
            else:
                logger.info("Report      : None")
        # Try to create a PDF if requested
        if config.make_pdf:
            _write_pdf()

    if config.export_plots:
        # Copy across the static plot images if requested
        _move_exported_plots()
        logger.info(
            "Plots       : {}{}".format(
                os.path.relpath(config.plots_dir),
                "   (overwritten)" if "export_plots" in overwritten else "",
            )
        )

    # Clean up temporary directory
    if clean_up:
        shutil.rmtree(report.tmp_dir)

    # Zip the data directory if requested
    if config.zip_data_dir and config.data_dir is not None:
        shutil.make_archive(config.data_dir, "zip", config.data_dir)
        shutil.rmtree(config.data_dir)


def _set_output_paths():
    """
    Set config output paths
    """

    # Add an output subdirectory if specified by template
    template_mod = config.avail_templates[config.template].load()
    try:
        config.output_dir = os.path.join(config.output_dir, template_mod.output_subdir)
    except AttributeError:
        pass  # No subdirectory variable given

    filename = config.filename

    if filename == "stdout":
        config.output_fn = sys.stdout
        logger.info("Printing report to stdout")
    else:
        if filename is not None and filename.endswith(".html"):
            filename = filename[:-5]
        if filename is None and config.title is not None:
            filename = re.sub(r"[^\w.-]", "", re.sub(r"[-\s]+", "-", config.title)).strip()
            filename += "_multiqc_report"
        if filename is not None:
            if "output_fn_name" not in config.nondefault_config:
                config.output_fn_name = f"{filename}.html"
            if "data_dir_name" not in config.nondefault_config:
                config.data_dir_name = f"{filename}_data"
            if "plots_dir_name" not in config.nondefault_config:
                config.plots_dir_name = f"{filename}_plots"
        if not config.output_fn_name.endswith(".html"):
            config.output_fn_name = f"{config.output_fn_name}.html"

        if config.make_report:
            config.output_fn = os.path.join(config.output_dir, config.output_fn_name)
        else:
            config.output_fn = None
        config.data_dir = os.path.join(config.output_dir, config.data_dir_name)
        config.plots_dir = os.path.join(config.output_dir, config.plots_dir_name)


def _move_exported_plots():
    """
    Assuming plots already exported to config.plots_tmp_dir(), move them to config.plots_dir
    """

    config.plots_dir = os.path.join(config.output_dir, config.plots_dir_name)
    if os.path.exists(config.plots_dir):
        if config.force:
            shutil.rmtree(config.plots_dir)
        else:
            logger.error(f"Output directory {config.plots_dir} already exists.")
            logger.info("Use -f or --force to overwrite existing reports")
            # if os.path.isdir(report.tmp_dir):
            #     shutil.rmtree(report.tmp_dir)
            raise RunError()

    # Modules have run, so plots directory should be complete by now. Move its contents.
    logger.debug(f"Moving plots directory from '{report.plots_tmp_dir()}' to '{config.plots_dir}'")
    shutil.copytree(
        report.plots_tmp_dir(),
        config.plots_dir,
        # Override default shutil.copy2 function to copy files. The default
        # function copies times and mode, which we want to avoid on purpose
        # to get around the problem with mounted CIFS shares (see #625).
        # shutil.copyfile only copies the file without any metadata.
        copy_function=shutil.copyfile,
    )
    shutil.rmtree(report.plots_tmp_dir())


def _create_or_override_dirs() -> Set[str]:
    """
    Write the report data to the output directory
    """
    overwritten = set()

    # Check for existing reports and remove if -f was specified
    if (
        (config.make_report and os.path.exists(config.output_fn))
        or (config.make_data_dir and os.path.exists(config.data_dir))
        or (config.export_plots and os.path.exists(config.plots_dir))
    ):
        if config.force:
            if config.make_report and os.path.exists(config.output_fn):
                overwritten.add("report")
                os.remove(config.output_fn)
            if config.make_data_dir and os.path.exists(config.data_dir):
                overwritten.add("data_dir")
                shutil.rmtree(config.data_dir)
            if config.export_plots and os.path.exists(config.plots_dir):
                overwritten.add("export_plots")
                shutil.rmtree(config.plots_dir)
        else:
            # Set up the base names of the report and the data dir
            report_num = 1
            dir_base = os.path.basename(config.data_dir)
            plots_base = os.path.basename(config.plots_dir)

            # Iterate through appended numbers until we find one that's free
            while (
                (config.make_report and os.path.exists(config.output_fn))
                or (config.make_data_dir and os.path.exists(config.data_dir))
                or (config.export_plots and os.path.exists(config.plots_dir))
            ):
                if config.make_report:
                    report_base, report_ext = os.path.splitext(config.output_fn_name)
                    config.output_fn = os.path.join(config.output_dir, f"{report_base}_{report_num}{report_ext}")
                config.data_dir = os.path.join(config.output_dir, f"{dir_base}_{report_num}")
                config.plots_dir = os.path.join(config.output_dir, f"{plots_base}_{report_num}")
                report_num += 1
            if config.make_report:
                config.output_fn_name = os.path.basename(config.output_fn)
            config.data_dir_name = os.path.basename(config.data_dir)
            config.plots_dir_name = os.path.basename(config.plots_dir)
            logger.info("Existing reports found, adding suffix to filenames. Use '--force' to overwrite.")

    if config.make_report:
        os.makedirs(os.path.dirname(config.output_fn), exist_ok=True)

    return overwritten


def _order_modules_and_sections():
    """
    Finalise modules and sections: place in the write order, add special-case modules
    """

    # In case if user passed exclude_modules or include_modules again:
    mod_anchors = include_or_exclude_modules([mod.anchor for mod in report.modules])
    for mod in report.modules:
        if mod.anchor not in mod_anchors:
            mod.hidden = True

    # Add section for software versions if any are found
    if not config.skip_versions_section and report.software_versions:
        # Importing here to avoid circular imports
        from multiqc.modules.software_versions import MultiqcModule as SoftwareVersionsModule

        # if the software versions module is not in report.modules, add it:
        if not any([isinstance(m, SoftwareVersionsModule) for m in report.modules]):
            report.modules.append(SoftwareVersionsModule())

    # Special-case module if we want to profile the MultiQC running time
    if config.profile_runtime:
        from multiqc.modules.profile_runtime import MultiqcModule as ProfileRuntimeModule

        # if the software versions module is not in report.modules, add it:
        if not any([isinstance(m, ProfileRuntimeModule) for m in report.modules]):
            report.modules.append(ProfileRuntimeModule())

    # Sort the report module output if we have a config
    if len(config.report_section_order) > 0:
        section_id_order = {}
        idx = 10
        for mod in reversed(report.modules):
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
        report.modules = [mod for i in reversed(sorted_ids) for mod in report.modules if mod.anchor == i]

    # Sort the report sections if we have a config
    # Basically the same as above, but sections within a module
    if len(config.report_section_order) > 0:
        # Go through each module
        for midx, mod in enumerate(report.modules):
            section_id_order = dict()
            # Get a list of the section anchors
            idx = 10
            for s in mod.sections:
                section_id_order[s.anchor] = idx
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
            report.modules[midx].sections = [s for i in sorted_ids for s in mod.sections if s.anchor == i]


def render_and_export_plots():
    """
    Render plot HTML, write PNG/SVG and plot data TSV/JSON to plots_tmp_dir() and data_tmp_dir(). Populates report.plot_data
    """
    for mod in report.modules:
        if mod.hidden:
            continue
        for s in mod.sections:
            if s.plot_id:
                plot = report.plot_by_id[s.plot_id]
                if isinstance(plot, Plot):
                    s.plot = plot.add_to_report(report)
                elif isinstance(plot, str):
                    s.plot = plot
                else:
                    logger.error(f"Unknown plot type for {mod.name} - {s.name}")
            else:
                s.plot = ""


def _render_general_stats_table() -> None:
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
            if "rid" not in h[k]:
                h[k]["rid"] = re.sub(r"\W+", "_", k).strip().strip("_")
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
        pconfig = {
            "id": "general_stats_table",
            "table_title": "General Statistics",
            "save_file": True,
            "raw_data_fn": "multiqc_general_stats",
        }
        p = table.plot(report.general_stats_data, report.general_stats_headers, pconfig)
        report.general_stats_html = p.add_to_report(clean_html_id=False)
    else:
        config.skip_generalstats = True


def _write_data_files() -> None:
    """
    Write auxiliary data files: module data exports, JSON dump, sources, dev plots data, upload MegaQC
    """
    # Modules have run, so data directory should be complete by now. Move its contents.
    logger.debug(f"Moving data file from '{report.data_tmp_dir()}' to '{config.data_dir}'")
    shutil.copytree(
        report.data_tmp_dir(),
        config.data_dir,
        # Override default shutil.copy2 function to copy files. The default
        # function copies times and mode, which we want to avoid on purpose
        # to get around the problem with mounted CIFS shares (see #625).
        # shutil.copyfile only copies the file without any metadata.
        copy_function=shutil.copyfile,
    )
    shutil.rmtree(report.data_tmp_dir())

    # Write the report sources to disk
    report.data_sources_tofile()

    # Create a file with the module DOIs
    report.dois_tofile(report.modules)

    # Data Export / MegaQC integration - save report data to file or send report data to an API endpoint
    if config.data_dump_file or (config.megaqc_url and config.megaqc_upload):
        dump = report.multiqc_dump_json()
        if config.data_dump_file:
            with (Path(config.data_dir) / "multiqc_data.json").open("w") as f:
                util_functions.dump_json(dump, f, indent=4, ensure_ascii=False)
        if config.megaqc_url:
            megaqc.multiqc_api_post(dump)

    if config.development:
        with open(os.path.join(config.data_dir, "multiqc_plots.js"), "w") as f:
            util_functions.dump_json(report.plot_data, f)


def _write_report():
    """
    Render and write report HTML to disk
    """
    # Copy over css & js files if requested by the theme
    for mod in report.modules:
        if mod.hidden:
            continue
        try:
            for to, path in mod.css.items():
                copy_to = os.path.join(report.tmp_dir, to)
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
            for to, path in mod.js.items():
                copy_to = os.path.join(report.tmp_dir, to)
                os.makedirs(os.path.dirname(copy_to))
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
        shutil.copytree(parent_template.template_dir, report.tmp_dir, dirs_exist_ok=True)

    # Copy the template files to the tmp directory (`dirs_exist_ok` makes sure
    # parent template files are overwritten)
    shutil.copytree(template_mod.template_dir, report.tmp_dir, dirs_exist_ok=True)

    # Function to include file contents in Jinja template
    def include_file(name, fdir=report.tmp_dir, b64=False):
        try:
            if fdir is None:
                fdir = ""
            path = os.path.join(fdir, name)

            if config.development:
                if os.path.exists(dev_path := os.path.join(template_mod.template_dir, name)):
                    fdir = template_mod.template_dir
                    name = dev_path
                    path = dev_path
                elif parent_template and os.path.exists(dev_path := os.path.join(parent_template.template_dir, name)):
                    fdir = template_mod.template_dir
                    name = dev_path
                    path = dev_path

                if re.match(r".*\.min\.(js|css)$", name):
                    unminimized_name = re.sub(r"\.min\.", ".", name)
                    if os.path.exists(os.path.join(fdir, unminimized_name)):
                        name = unminimized_name

                if name.endswith(".js"):
                    return f'</script><script type="text/javascript" src="{name}">'
                if name.endswith(".css"):
                    return f'</style><link rel="stylesheet" href="{name}">'

            if b64:
                with io.open(path, "rb") as f:
                    return base64.b64encode(f.read()).decode("utf-8")
            else:
                with io.open(path, "r", encoding="utf-8") as f:
                    return f.read()
        except (OSError, IOError) as e:
            logger.error(f"Could not include file '{name}': {e}")

    # Load the report template
    try:
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(report.tmp_dir))
        env.globals["include_file"] = include_file
        j_template = env.get_template(template_mod.base_fn, globals={"development": config.development})
    except:  # noqa: E722
        raise IOError(f"Could not load {config.template} template file '{template_mod.base_fn}'")

    # Compress the report plot JSON data
    runtime_compression_start = time.time()
    logger.debug("Compressing plot data")
    report.plot_compressed_json = report.compress_json(report.plot_data)
    report.runtimes["total_compression"] = time.time() - runtime_compression_start

    # Use jinja2 to render the template and overwrite
    report.analysis_files = [os.path.realpath(d) for d in report.analysis_files]
    report_output = j_template.render(report=report, config=config)
    if config.filename == "stdout":
        print(report_output.encode("utf-8"), file=sys.stdout)
    else:
        try:
            with io.open(config.output_fn, "w", encoding="utf-8") as f:
                print(report_output, file=f)
        except IOError as e:
            raise IOError(f"Could not print report to '{config.output_fn}' - {IOError(e)}")

        # Copy over files if requested by the theme
        try:
            for f in template_mod.copy_files:
                fn = os.path.join(report.tmp_dir, f)
                dest_dir = os.path.join(os.path.dirname(config.output_fn), f)
                shutil.copytree(fn, dest_dir, dirs_exist_ok=True)
        except AttributeError:
            pass  # No files to copy


def _write_pdf():
    try:
        pdf_fn_name = config.output_fn.replace(".html", ".pdf")
        pandoc_call = [
            "pandoc",
            "--standalone",
            config.output_fn,
            "--output",
            pdf_fn_name,
            "--pdf-engine=xelatex",
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
        pdf_exit_code = subprocess.call(pandoc_call)
        if pdf_exit_code != 0:
            logger.error("Error creating PDF! Pandoc returned a non-zero exit code.")
        else:
            logger.info(f"PDF Report  : {pdf_fn_name}")
    except OSError as e:
        if e.errno == errno.ENOENT:
            logger.error("Error creating PDF - pandoc not found. Is it installed? http://pandoc.org/")
        else:
            logger.error(
                "Error creating PDF! Something went wrong when creating the PDF\n"
                + ("=" * 60)
                + f"\n{traceback.format_exc()}\n"
                + ("=" * 60)
            )
