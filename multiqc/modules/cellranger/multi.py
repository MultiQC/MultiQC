"""MultiQC module to parse output from Cell Ranger multi"""

import json
import logging
import re
from typing import Dict, Tuple, List

from multiqc import BaseMultiqcModule
from multiqc.plots import table, linegraph
from multiqc.modules.cellranger.utils import (
    resolve_dict,
    set_hidden_cols,
    SectionData,
    update_data_and_header,
    combine_data_and_headers,
)

# Initialise the logger
log = logging.getLogger(__name__)


def parse_multi_html(module: BaseMultiqcModule):
    """Cell Ranger multi report parser"""

    def add_gex_sections():
        if len(gex_alerts.data) > 0:
            module.add_section(
                name="Gene-Expression alerts",
                anchor="cellranger-multi-gex-alerts",
                plot=table.plot(
                    data=gex_alerts.data,
                    headers=gex_alerts.headers,
                    pconfig={
                        "namespace": "Multi",
                        "id": "cellranger-multi-gex-alerts-table",
                        "title": "CellRanger Multi: Gene-Expression alerts",
                    },
                ),
            )

        if len(gex_cells_metrics.data) > 0:
            module.add_section(
                name="Gene-Expression cell metrics",
                anchor="cellranger-multi-gex-cell-metrics",
                plot=table.plot(
                    data=gex_cells_metrics.data,
                    headers=gex_cells_metrics.headers,
                    pconfig={
                        "namespace": "Multi",
                        "id": "cellranger-multi-gex-cell-table",
                        "title": "CellRanger Multi: Gene-Expression cell metrics",
                    },
                ),
            )

        if len(gex_library_metrics_summary.data) > 0:
            module.add_section(
                name="Gene-Expression Library Metrics Summary",
                anchor="cellranger-multi-gex-library-metrics",
                plot=table.plot(
                    data=gex_library_metrics_summary.data,
                    headers=set_hidden_cols(
                        gex_library_metrics_summary.headers,
                        [
                            "Number of short reads skipped",
                            "Number of reads",
                            "Fastq ID",
                            "Physical library ID",
                            "Mean reads per cell",
                        ],
                    ),
                    pconfig={
                        "namespace": "Multi",
                        "id": "cellranger-multi-gex-library-table",
                        "title": "CellRanger Multi: Gene-Expression library metrics summary",
                    },
                ),
            )

        if len(gex_bc_plot.data) > 0:
            module.add_section(
                name="Gene-Expression Barcode Rank Plot",
                anchor="cellranger-multi-gex-bcrank-section",
                helptext=gex_bc_plot.headers["helptext"],
                plot=linegraph.plot(gex_bc_plot.data, gex_bc_plot.headers["config"]),
            )

        if len(gex_genes_plot.data) > 0:
            module.add_section(
                name="Gene-Expression Median Genes per Cell",
                anchor="cellranger-multi-gex-genes-section",
                helptext=gex_genes_plot.headers["helptext"],
                plot=linegraph.plot(gex_genes_plot.data, gex_genes_plot.headers["config"]),
            )

    def add_vdj_t_sections():
        if len(vdj_t_alerts.data) > 0:
            module.add_section(
                name="VDJ-T alerts",
                anchor="cellranger-multi-vdj-t-alerts",
                plot=table.plot(
                    data=vdj_t_alerts.data,
                    headers=vdj_t_alerts.headers,
                    pconfig={
                        "namespace": "Multi",
                        "id": "cellranger-multi-vdj-t-alerts-table",
                        "title": "CellRanger Multi: VDJ-T alerts",
                    },
                ),
            )

        if len(vdj_t_expression_annotation_metrics.data) > 0:
            module.add_section(
                name="VDJ-T expression and annotation metrics",
                anchor="cellranger-multi-vdj-t-expression-annotation-metrics",
                plot=table.plot(
                    data=vdj_t_expression_annotation_metrics.data,
                    headers=set_hidden_cols(
                        vdj_t_expression_annotation_metrics.headers,
                        [
                            "Median TRA UMIs per Cell",
                            "Median TRB UMIs per Cell",
                            "Number of cells with productive V-J spanning pair",
                        ],
                    ),
                    pconfig={
                        "namespace": "Multi",
                        "id": "cellranger-multi-vdj-t-expression-annotation-table",
                        "title": "CellRanger Multi: VDJ-T alerts",
                    },
                ),
            )

        if len(vdj_t_bc_plot.data) > 0:
            module.add_section(
                name="VDJ-T Barcode Rank Plot",
                anchor="cellranger-multi-vdj-t-bcrank-section",
                helptext=vdj_t_bc_plot.headers["helptext"],
                plot=linegraph.plot(vdj_t_bc_plot.data, vdj_t_bc_plot.headers["config"]),
            )

    def add_vdj_b_sections():
        if len(vdj_b_alerts.data) > 0:
            module.add_section(
                name="VDJ-B alerts",
                anchor="cellranger-multi-vdj-b-alerts",
                plot=table.plot(
                    data=vdj_b_alerts.data,
                    headers=vdj_b_alerts.headers,
                    pconfig={
                        "namespace": "Multi",
                        "id": "cellranger-multi-vdj-b-alerts-table",
                        "title": "CellRanger Multi: VDJ-B alerts",
                    },
                ),
            )

        if len(vdj_b_expression_annotation_metrics.data) > 0:
            module.add_section(
                name="VDJ-B expression and annotation metrics",
                anchor="cellranger-multi-vdj-b-expression-annotation-metrics",
                plot=table.plot(
                    data=vdj_b_expression_annotation_metrics.data,
                    headers=set_hidden_cols(
                        vdj_b_expression_annotation_metrics.headers,
                        [
                            "Median TRA UMIs per Cell",
                            "Median TRB UMIs per Cell",
                            "Number of cells with productive V-J spanning pair",
                        ],
                    ),
                    pconfig={
                        "namespace": "Multi",
                        "id": "cellranger-multi-vdj-b-expression-annotation-table",
                        "title": "CellRanger Multi: VDJ-B alerts",
                    },
                ),
            )

        if len(vdj_b_bc_plot.data) > 0:
            module.add_section(
                name="VDJ-B Barcode Rank Plot",
                anchor="cellranger-multi-vdj-b-bcrank-section",
                helptext=vdj_b_bc_plot.headers["helptext"],
                plot=linegraph.plot(vdj_b_bc_plot.data, vdj_b_bc_plot.headers["config"]),
            )

    def add_antibody_sections():
        if len(antibody_alerts.data) > 0:
            module.add_section(
                name="Antibody alerts",
                anchor="cellranger-multi-antibody-alerts",
                plot=table.plot(
                    data=antibody_alerts.data,
                    headers=antibody_alerts.headers,
                    pconfig={
                        "namespace": "Multi",
                        "id": "cellranger-multi-antibody-alerts-table",
                        "title": "CellRanger Multi: Antibody alerts",
                    },
                ),
            )

        if len(antibody_expression_mapping_metrics.data) > 0:
            module.add_section(
                name="Antibody expression and mapping metrics",
                anchor="cellranger-multi-antibody-expression-mapping-metrics",
                plot=table.plot(
                    data=antibody_expression_mapping_metrics.data,
                    headers=set_hidden_cols(
                        antibody_expression_mapping_metrics.headers, ["Mean antibody reads usable per cell"]
                    ),
                    pconfig={
                        "namespace": "Multi",
                        "id": "cellranger-multi-antibody-expression-mapping-table",
                        "title": "CellRanger Multi: VDJ-B alerts",
                    },
                ),
            )

        if len(antibody_bc_plot.data) > 0:
            module.add_section(
                name="Antibody Barcode Rank Plot",
                anchor="cellranger-multi-antibody-bcrank-section",
                helptext=antibody_bc_plot.headers["helptext"],
                plot=linegraph.plot(antibody_bc_plot.data, antibody_bc_plot.headers["config"]),
            )

    general_data = SectionData()
    gex_alerts = SectionData()
    gex_cells_metrics = SectionData()
    gex_library_metrics_summary = SectionData()
    gex_bc_plot = SectionData()
    gex_genes_plot = SectionData()
    vdj_t_alerts = SectionData()
    vdj_t_expression_annotation_metrics = SectionData()
    vdj_t_bc_plot = SectionData()
    vdj_b_alerts = SectionData()
    vdj_b_expression_annotation_metrics = SectionData()
    vdj_b_bc_plot = SectionData()
    antibody_alerts = SectionData()
    antibody_expression_mapping_metrics = SectionData()
    antibody_bc_plot = SectionData()

    gex_bc_plot.headers = {
        "helptext": "",
        "config": {
            "id": "cellranger-multi-gex-bcr-plot",
            "title": "Barcode Rank Plot",
            "xlab": "Barcodes",
            "ylab": "UMI counts",
            "xlog": True,
            "ylog": True,
            "xmin": 1,
            "ymin": 1,
        },
    }
    gex_genes_plot.headers = {
        "helptext": "",
        "config": {
            "id": "cellranger-multi-gex-genes-plot",
            "title": "Median Genes per Cell",
            "xlab": "Mean Reads per Cell",
            "ylab": "Median Genes per Cell",
        },
    }
    vdj_t_bc_plot.headers = {
        "helptext": "",
        "config": {
            "id": "cellranger-multi-vdj-t-bcr-plot",
            "title": "V(D)J Barcode Rank Plot",
            "xlab": "Barcodes",
            "ylab": "UMI counts",
            "xlog": True,
            "ylog": True,
            "xmin": 1,
            "ymin": 1,
        },
    }
    vdj_b_bc_plot.headers = {
        "helptext": "",
        "config": {
            "id": "cellranger-multi-vdj-b-bcr-plot",
            "title": "V(D)J Barcode Rank Plot",
            "xlab": "Barcodes",
            "ylab": "UMI counts",
            "xlog": True,
            "ylog": True,
            "xmin": 1,
            "ymin": 1,
        },
    }
    antibody_bc_plot.headers = {
        "helptext": "",
        "config": {
            "id": "cellranger-multi-antibody-bcr-plot",
            "title": "AB Barcode Rank Plot",
            "xlab": "Barcodes",
            "ylab": "UMI counts",
            "xlog": True,
            "ylog": True,
            "xmin": 1,
            "ymin": 1,
        },
    }

    for file in module.find_log_files("cellranger/multi_html", filehandles=True):
        data = {}
        for line in file["f"]:
            line = line.strip()
            if line.startswith("const data"):
                line = line.replace("const data = ", "")
                data = json.loads(line)
                break

        sample_name = module.clean_s_name(data["sample"]["id"], file)

        try:
            version = data["data"]["sample_websummary"]["header_info"]["Pipeline Version"]
            version_match = re.search(r"cellranger-([\d\.]+)", version)
            if version_match:
                module.add_software_version(version_match.group(1), sample_name)
        except (KeyError, AssertionError):
            log.debug(f"Unable to parse version for sample {sample_name}")

        _build_gex_data(
            data["data"],
            sample_name,
            gex_alerts,
            gex_cells_metrics,
            gex_library_metrics_summary,
            gex_bc_plot,
            gex_genes_plot,
        )

        _build_vdj_t_data(data["data"], sample_name, vdj_t_alerts, vdj_t_expression_annotation_metrics, vdj_t_bc_plot)

        _build_vdj_b_data(data["data"], sample_name, vdj_b_alerts, vdj_b_expression_annotation_metrics, vdj_b_bc_plot)

        _build_antibody_data(
            data["data"], sample_name, antibody_alerts, antibody_expression_mapping_metrics, antibody_bc_plot
        )

        module.add_data_source(file, sample_name, module="cellranger", section="multi")

    general_data.data, general_data.headers = combine_data_and_headers(
        [
            (gex_cells_metrics.data, gex_cells_metrics.headers, ["Cells", "Mean reads per Cell"]),
            (
                gex_library_metrics_summary.data,
                gex_library_metrics_summary.headers,
                ["Q30 barcodes", "Sequencing saturation"],
            ),
        ]
    )

    general_data.data = module.ignore_samples(general_data.data)
    if len(general_data.data) == 0:
        return 0
    module.general_stats_addcols(general_data.data, general_data.headers)

    module.write_data_file(gex_cells_metrics.data, "multiqc_cellranger_gex_cells")
    module.write_data_file(gex_library_metrics_summary.data, "multiqc_cellranger_gex_libraries")
    module.write_data_file(vdj_t_expression_annotation_metrics.data, "multiqc_cellranger_vdj_t_expr_anno")
    module.write_data_file(vdj_b_expression_annotation_metrics.data, "multiqc_cellranger_vdj_b_expr_anno")
    module.write_data_file(antibody_expression_mapping_metrics.data, "multiqc_cellranger_antibody_expr_mapp")

    for section in [
        gex_alerts,
        gex_cells_metrics,
        gex_library_metrics_summary,
        gex_bc_plot,
        gex_genes_plot,
        vdj_t_alerts,
        vdj_t_expression_annotation_metrics,
        vdj_t_bc_plot,
        vdj_b_alerts,
        vdj_b_expression_annotation_metrics,
        vdj_b_bc_plot,
        antibody_alerts,
        antibody_expression_mapping_metrics,
        antibody_bc_plot,
    ]:
        section.data = module.ignore_samples(section.data)

    add_gex_sections()
    add_vdj_t_sections()
    add_vdj_b_sections()
    add_antibody_sections()
    return len(general_data.data)


def _parse_table(table_root: Dict, color_dict=None) -> Tuple[Dict, Dict]:
    if not table_root:
        return {}, {}

    if color_dict is None:
        color_dict = {}

    if len(table_root["table"]["rows"]) > 1:
        log.critical(f'Multiple rows for table: {table_root["help"].get("title", "UNDEFINED")}')

    table_keys = table_root["table"]["header"]
    table_data = table_root["table"]["rows"][0]
    table_help = table_root["help"]

    parsed_data = dict(zip(table_keys, table_data))
    parsed_headers = {entry[0]: {"title": entry[0], "description": entry[1][0]} for entry in table_help["data"]}

    parsed_data, parsed_headers = update_data_and_header(parsed_data, parsed_headers, color_dict)
    return parsed_data, parsed_headers


def _parse_alerts(data: Dict, path: List):
    try:
        alert_list = resolve_dict(data, path)
    except KeyError:
        alert_list = []
    alerts = {alert["title"]: alert["level"] if alert["level"] != "ERROR" else "FAIL" for alert in alert_list}
    return alerts


def _build_gex_data(
    data: Dict,
    sample_name: str,
    alerts: SectionData,
    cells_metrics: SectionData,
    library_metrics: SectionData,
    barcode_plot: SectionData,
    genes_plot: SectionData,
):
    sample_websummary = data.get("sample_websummary", {}).get("gex_tab", {})
    library_websummary = data.get("library_websummary", {}).get("gex_tab", {})

    if not sample_websummary or not library_websummary:
        log.debug("No Gene-Expression data found")
        return
    if not sample_websummary and not library_websummary:
        log.debug("Gene-Expression data incomplete")
        return

    parsed_alerts = _parse_alerts(sample_websummary, ["alerts"])
    parsed_alerts.update(_parse_alerts(library_websummary, ["alerts"]))
    if parsed_alerts:
        alerts.update_sample(parsed_alerts, {}, sample_name)

    try:
        cells_source = sample_websummary["content"]["hero_metrics"]
        cells, cells_headers = _parse_table(
            cells_source,
            {
                "Cells": "GnBu",
                "Mean reads per cell": "BuPu",
                "Median genes per Cell": "Purples",
                "Total genes detected": "GnBu",
                "Median UMI counts per Cell": "BuPu",
                "Confidently mapped reads in cells": "Purples",
            },
        )

        if cells:
            cells_metrics.update_sample(cells, cells_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/gex_tab/content/hero_metrics")

    try:
        sequencing_metrics_source = library_websummary["content"]["sequencing_metrics_table"]
        sequencing_metrics, sequencing_metrics_headers = _parse_table(
            sequencing_metrics_source,
            {"Q30 barcodes": "GnBu", "Q30 UMI": "BuPu", "Q30 RNA read": "Purples", "Q30 RNA read2": "GnBu"},
        )

        if sequencing_metrics:
            library_metrics.update_sample(sequencing_metrics, sequencing_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/gex_tab/content/sequencing_metrics_table")

    try:
        physical_library_metrics_source = library_websummary["content"]["physical_library_metrics_table"]
        physical_library_metrics, physical_library_metrics_headers = _parse_table(
            physical_library_metrics_source,
            {"Valid barcodes": "GnBu", "Valid UMIs": "BuPu", "Sequencing saturation": "Purples"},
        )

        if physical_library_metrics:
            library_metrics.update_sample(physical_library_metrics, physical_library_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/gex_tab/content/physical_library_metrics_table")

    try:
        bc_plot_source = library_websummary["content"]["barcode_rank_plot"]
        bc_plot_data = {}
        for subset in bc_plot_source["plot"]["data"]:
            bc_plot_data.update({x: y for x, y in zip(subset["x"], subset["y"])})

        if bc_plot_data:
            barcode_plot.update_sample(bc_plot_data, {"helptext": bc_plot_source["help"]["helpText"]}, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/gex_tab/content/barcode_rank_plot")

    try:
        genes_plot_source = library_websummary["content"]["median_genes_per_cell_plot"]
        genes_plot_data = {}
        for subset in genes_plot_source["plot"]["data"]:
            genes_plot_data.update({x: y for x, y in zip(subset["x"], subset["y"])})

        if genes_plot_data:
            genes_plot.update_sample(genes_plot_data, {"helptext": genes_plot_source["help"]["helpText"]}, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/gex_tab/content/median_genes_per_cell_plot")


def _build_vdj_t_data(
    data: Dict,
    sample_name: str,
    alerts: SectionData,
    expression_annotation_metrics: SectionData,
    barcode_plot: SectionData,
):
    sample_websummary = data.get("sample_websummary", {}).get("vdj_t_tab", {})
    library_websummary = data.get("library_websummary", {}).get("vdj_t_tab", {})

    if not sample_websummary and not library_websummary:
        log.debug("No VDJ-T data found")
        return
    if not sample_websummary or not library_websummary:
        log.debug("VDJ-T data incomplete")
        return

    parsed_alerts = _parse_alerts(sample_websummary, ["alerts"])
    parsed_alerts.update(_parse_alerts(library_websummary, ["alerts"]))
    if parsed_alerts:
        alerts.update_sample(parsed_alerts, {}, sample_name)

    try:
        expression_metrics_source = sample_websummary["content"]["hero_metrics"]
        expression_metrics, expression_metrics_headers = _parse_table(
            expression_metrics_source,
            {
                "Estimated number of cells": "GnBu",
                "Number of cells with productive V-J spanning pair": "BuPu",
                "Median TRA UMIs per Cell": "Purples",
                "Median TRB UMIs per Cell": "GnBu",
            },
        )

        if expression_metrics:
            expression_annotation_metrics.update_sample(expression_metrics, expression_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/vdj_t_tab/content/hero_metrics")

    try:
        annotation_metrics_source = sample_websummary["content"]["annotation_metrics_table"]
        annotation_metrics, annotation_metrics_headers = _parse_table(
            annotation_metrics_source,
            {
                "Cells with productive V-J spanning pair": "GnBu",
                "Cells with productive V-J spanning (TRA, TRB) pair": "BuPu",
                "Cells with productive TRA contig": "Purples",
                "Cells with productive TRB contig": "GnBu",
                "Paired clonotype diversity": "BuPu",
            },
        )

        if annotation_metrics:
            expression_annotation_metrics.update_sample(annotation_metrics, annotation_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/vdj_t_tab/content/annotation_metrics_table")

    try:
        bc_plot_source = library_websummary["content"]["barcode_rank_plot"]
        bc_plot_data = {}
        for subset in bc_plot_source["plot"]["data"]:
            bc_plot_data.update(dict(zip(subset["x"], subset["y"])))

        if bc_plot_data:
            barcode_plot.update_sample(bc_plot_data, {"helptext": bc_plot_source["help"]["helpText"]}, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/vdj_t_tab/content/barcode_rank_plot")


def _build_vdj_b_data(
    data: Dict,
    sample_name: str,
    alerts: SectionData,
    expression_annotation_metrics: SectionData,
    barcode_plot: SectionData,
):
    sample_websummary = data.get("sample_websummary", {}).get("vdj_b_tab", {})
    library_websummary = data.get("library_websummary", {}).get("vdj_b_tab", {})

    if not sample_websummary and not library_websummary:
        log.debug("No VDJ-B data found")
        return
    if not sample_websummary or not library_websummary:
        log.debug("VDJ-B data incomplete")
        return

    parsed_alerts = _parse_alerts(sample_websummary, ["alerts"])
    parsed_alerts.update(_parse_alerts(data, ["alerts"]))
    if parsed_alerts:
        alerts.update_sample(parsed_alerts, {}, sample_name)

    try:
        expression_metrics_source = sample_websummary["content"]["hero_metrics"]
        expression_metrics, expression_metrics_headers = _parse_table(
            expression_metrics_source,
            {
                "Estimated number of cells": "GnBu",
                "Number of cells with productive V-J spanning pair": "BuPu",
                "Median IGH UMIs per Cell": "Purples",
                "Median IGK UMIs per Cell": "GnBu",
                "Median IGL UMIs per Cell": "BuPu",
            },
        )

        if expression_metrics:
            expression_annotation_metrics.update_sample(expression_metrics, expression_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/vdj_b_tab/content/hero_metrics")

    try:
        annotation_metrics_source = sample_websummary["content"]["annotation_metrics_table"]
        annotation_metrics, annotation_metrics_headers = _parse_table(
            annotation_metrics_source,
            {
                "Cells with productive V-J spanning pair": "GnBu",
                "Cells with productive V-J spanning (IGK, IGH) pair": "BuPu",
                "Cells with productive V-J spanning (IGL, IGH) pair": "Purples",
                "Cells with productive IGH contig": "GnBu",
                "Cells with productive IGK contig": "BuPu",
                "Cells with productive IGL contig": "Purples",
                "Paired clonotype diversity": "GnBu",
            },
        )

        if annotation_metrics:
            expression_annotation_metrics.update_sample(annotation_metrics, annotation_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/vdj_b_tab/content/annotation_metrics_table")

    try:
        bc_plot_source = library_websummary["content"]["barcode_rank_plot"]
        bc_plot_data = {}
        for subset in bc_plot_source["plot"]["data"]:
            bc_plot_data.update(dict(zip(subset["x"], subset["y"])))

        if bc_plot_data:
            barcode_plot.update_sample(bc_plot_data, {"helptext": bc_plot_source["help"]["helpText"]}, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/vdj_b_tab/content/barcode_rank_plot")


def _build_antibody_data(
    data: Dict,
    sample_name: str,
    alerts: SectionData,
    expression_mapping_metrics: SectionData,
    barcode_plot: SectionData,
):
    sample_websummary = data.get("sample_websummary", {}).get("antibody_tab", {})
    library_websummary = data.get("library_websummary", {}).get("antibody_tab", {})

    if not sample_websummary and not library_websummary:
        log.debug("No antibody data found")
        return ()
    if not sample_websummary or not library_websummary:
        log.debug("Antibody data incomplete")
        return ()

    parsed_alerts = _parse_alerts(sample_websummary, ["alerts"])
    parsed_alerts.update(_parse_alerts(data, ["alerts"]))
    if parsed_alerts:
        alerts.update_sample(parsed_alerts, {}, sample_name)

    try:
        expression_metrics_source = sample_websummary["content"]["hero_metrics"]
        expression_metrics, expression_metrics_headers = _parse_table(
            expression_metrics_source,
            {
                "Cells": "GnBu",
                "Median UMI counts per cell": "BuPu",
                "Mean antibody reads usable per cell": "Purples",
                "Antibody reads in cells": "GnBu",
            },
        )

        if expression_metrics:
            expression_mapping_metrics.update_sample(expression_metrics, expression_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/antibody_tab/content/hero_metrics")

    try:
        bc_plot_source = library_websummary["content"]["barcode_rank_plot"]
        bc_plot_data = {}
        for subset in bc_plot_source["plot"]["data"]:
            bc_plot_data.update(dict(zip(subset["x"], subset["y"])))

        if bc_plot_data:
            barcode_plot.update_sample(bc_plot_data, {"helptext": bc_plot_source["help"]["helpText"]}, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/antibody_tab/content/barcode_rank_plot")
