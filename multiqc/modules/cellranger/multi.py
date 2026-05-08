"""MultiQC module to parse output from Cell Ranger multi"""

import json
import logging
import re
from typing import Dict, Tuple, List

from multiqc import BaseMultiqcModule
from multiqc.plots import table, linegraph
from multiqc.modules.cellranger.utils import (
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
                name="Gene-Expression: Alerts",
                anchor="cellranger-multi-gex-alerts",
                plot=table.plot(
                    data=gex_alerts.data,
                    headers=set_hidden_cols(gex_alerts.headers, ["Intron mode used"]),
                    pconfig={
                        "namespace": "Multi-Gex",
                        "id": "cellranger-multi-gex-alerts-table",
                        "title": "CellRanger Multi: Gene-Expression Alerts",
                    },
                ),
            )

        if len(gex_cells_metrics.data) > 0:
            module.add_section(
                name="Gene-Expression: Cell Metrics",
                anchor="cellranger-multi-gex-cell-metrics",
                plot=table.plot(
                    data=gex_cells_metrics.data,
                    headers=set_hidden_cols(gex_cells_metrics.headers, ["Median UMI counts per cell"]),
                    pconfig={
                        "namespace": "Multi-Gex",
                        "id": "cellranger-multi-gex-cell-table",
                        "title": "CellRanger Multi: Gene-Expression Cell Metrics",
                    },
                ),
            )

        if len(gex_library_metrics_summary.data) > 0:
            module.add_section(
                name="Gene-Expression: Library Metrics Summary",
                anchor="cellranger-multi-gex-library-metrics",
                plot=table.plot(
                    data=gex_library_metrics_summary.data,
                    headers=set_hidden_cols(
                        gex_library_metrics_summary.headers,
                        [
                            "Confidently mapped reads in cells",
                            "Fastq ID",
                            "Number of reads",
                            "Number of short reads skipped",
                            "Physical library ID",
                            "Mean reads per cell",
                        ],
                    ),
                    pconfig={
                        "namespace": "Multi-Gex",
                        "id": "cellranger-multi-gex-library-table",
                        "title": "CellRanger Multi: Gene-Expression Library Metrics Summary",
                    },
                ),
            )

        if len(gex_bc_plot.data) > 0:
            module.add_section(
                name="Gene-Expression: Barcode Rank Plot",
                anchor="cellranger-multi-gex-bcrank-section",
                helptext=gex_bc_plot.headers["helptext"],
                plot=linegraph.plot(gex_bc_plot.data, gex_bc_plot.headers["config"]),
            )

        if len(gex_genes_plot.data) > 0:
            module.add_section(
                name="Gene-Expression: Median Genes per Cell",
                anchor="cellranger-multi-gex-genes-section",
                helptext=gex_genes_plot.headers["helptext"],
                plot=linegraph.plot(gex_genes_plot.data, gex_genes_plot.headers["config"]),
            )

        if len(gex_sequencing_plot.data) > 0:
            module.add_section(
                name="Gene-Expression: Sequencing Saturation",
                anchor="cellranger-multi-gex-sequencing-section",
                helptext=gex_sequencing_plot.headers["helptext"],
                plot=linegraph.plot(gex_sequencing_plot.data, gex_sequencing_plot.headers["config"]),
            )

    def add_vdj_t_sections():
        if len(vdj_t_alerts.data) > 0:
            module.add_section(
                name="VDJ-T: Alerts",
                anchor="cellranger-multi-vdj-t-alerts",
                plot=table.plot(
                    data=vdj_t_alerts.data,
                    headers=vdj_t_alerts.headers,
                    pconfig={
                        "namespace": "Multi-VDJ-T",
                        "id": "cellranger-multi-vdj-t-alerts-table",
                        "title": "CellRanger Multi: VDJ-T Alerts",
                    },
                ),
            )

        if len(vdj_t_expression_metrics.data) > 0:
            module.add_section(
                name="VDJ-T: Expression Metrics",
                anchor="cellranger-multi-vdj-t-expression-metrics",
                plot=table.plot(
                    data=vdj_t_expression_metrics.data,
                    headers=vdj_t_expression_metrics.headers,
                    pconfig={
                        "namespace": "Multi-VDJ-T",
                        "id": "cellranger-multi-vdj-t-expression-table",
                        "title": "CellRanger Multi: VDJ-T Expression Metrics",
                    },
                ),
            )

        if len(vdj_t_annotation_metrics.data) > 0:
            module.add_section(
                name="VDJ-T: Annotation Metrics",
                anchor="cellranger-multi-vdj-t-annotation-metrics",
                plot=table.plot(
                    data=vdj_t_annotation_metrics.data,
                    headers=set_hidden_cols(
                        vdj_t_annotation_metrics.headers,
                        ["Cells with productive V-J spanning (TRA, TRB) pair"],
                    ),
                    pconfig={
                        "namespace": "Multi-VDJ-T",
                        "id": "cellranger-multi-vdj-t-annotation-table",
                        "title": "CellRanger Multi: VDJ-T Annotation Metrics",
                    },
                ),
            )

        if len(vdj_t_bc_plot.data) > 0:
            module.add_section(
                name="VDJ-T: Barcode Rank Plot",
                anchor="cellranger-multi-vdj-t-bcrank-section",
                helptext=vdj_t_bc_plot.headers["helptext"],
                plot=linegraph.plot(vdj_t_bc_plot.data, vdj_t_bc_plot.headers["config"]),
            )

    def add_vdj_b_sections():
        if len(vdj_b_alerts.data) > 0:
            module.add_section(
                name="VDJ-B: Alerts",
                anchor="cellranger-multi-vdj-b-alerts",
                plot=table.plot(
                    data=vdj_b_alerts.data,
                    headers=vdj_b_alerts.headers,
                    pconfig={
                        "namespace": "Multi-VDJ-B",
                        "id": "cellranger-multi-vdj-b-alerts-table",
                        "title": "CellRanger Multi: VDJ-B Alerts",
                    },
                ),
            )

        if len(vdj_b_expression_metrics.data) > 0:
            module.add_section(
                name="VDJ-B: Expression Metrics",
                anchor="cellranger-multi-vdj-b-expression-metrics",
                plot=table.plot(
                    data=vdj_b_expression_metrics.data,
                    headers=vdj_b_expression_metrics.headers,
                    pconfig={
                        "namespace": "Multi-VDJ-B",
                        "id": "cellranger-multi-vdj-b-expression-table",
                        "title": "CellRanger Multi: VDJ-B Expression Metrics",
                    },
                ),
            )

        if len(vdj_b_annotation_metrics.data) > 0:
            module.add_section(
                name="VDJ-B: Annotation Metrics",
                anchor="cellranger-multi-vdj-b-annotation-metrics",
                plot=table.plot(
                    data=vdj_b_annotation_metrics.data,
                    headers=set_hidden_cols(
                        vdj_b_annotation_metrics.headers,
                        [
                            "Cells with productive V-J spanning (IGK, IGH) pair",
                            "Cells with productive V-J spanning (IGL, IGH) pair",
                        ],
                    ),
                    pconfig={
                        "namespace": "Multi-VDJ-B",
                        "id": "cellranger-multi-vdj-b-annotation-table",
                        "title": "CellRanger Multi: VDJ-B Annotation Metrics",
                    },
                ),
            )

        if len(vdj_b_bc_plot.data) > 0:
            module.add_section(
                name="VDJ-B: Barcode Rank Plot",
                anchor="cellranger-multi-vdj-b-bcrank-section",
                helptext=vdj_b_bc_plot.headers["helptext"],
                plot=linegraph.plot(vdj_b_bc_plot.data, vdj_b_bc_plot.headers["config"]),
            )

    def add_antibody_sections():
        if len(antibody_alerts.data) > 0:
            module.add_section(
                name="Antibody: alerts",
                anchor="cellranger-multi-antibody-alerts",
                plot=table.plot(
                    data=antibody_alerts.data,
                    headers=antibody_alerts.headers,
                    pconfig={
                        "namespace": "Multi-Antibody",
                        "id": "cellranger-multi-antibody-alerts-table",
                        "title": "CellRanger Multi: Antibody alerts",
                    },
                ),
            )

        if len(antibody_expression_mapping_metrics.data) > 0:
            module.add_section(
                name="Antibody: expression and mapping metrics",
                anchor="cellranger-multi-antibody-expression-mapping-metrics",
                plot=table.plot(
                    data=antibody_expression_mapping_metrics.data,
                    headers=set_hidden_cols(antibody_expression_mapping_metrics.headers, []),
                    pconfig={
                        "namespace": "Multi-Antibody",
                        "id": "cellranger-multi-antibody-expression-mapping-table",
                        "title": "CellRanger Multi: Antibody Expression and Mapping metrics",
                    },
                ),
            )

        if len(antibody_bc_plot.data) > 0:
            module.add_section(
                name="Antibody: Barcode Rank Plot",
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
    gex_sequencing_plot = SectionData()
    vdj_t_alerts = SectionData()
    vdj_t_expression_metrics = SectionData()
    vdj_t_annotation_metrics = SectionData()
    vdj_t_bc_plot = SectionData()
    vdj_b_alerts = SectionData()
    vdj_b_expression_metrics = SectionData()
    vdj_b_annotation_metrics = SectionData()
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
            "xmin": 0,
            "ymin": 0,
        },
    }
    gex_sequencing_plot.headers = {
        "helptext": "",
        "config": {
            "id": "cellranger-multi-gex-sequencing-plot",
            "title": "Sequencing Saturation",
            "xlab": "Mean Reads per Cell",
            "ylab": "Sequencing Saturation",
            "xmin": 0,
            "ymin": 0,
            "ymax": 1,
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

        crm_version = "UNDEFINED"
        try:
            version = data["data"]["sample_websummary"]["header_info"]["Pipeline Version"]
            version_match = re.search(r"cellranger-([\d\.]+)", version)
            if version_match:
                crm_version = version_match.group(1)
                module.add_software_version(crm_version, sample_name)
        except KeyError:
            log.debug(f"Unable to parse version for sample {sample_name}")

        _build_gex_data(
            data["data"],
            sample_name,
            gex_alerts,
            gex_cells_metrics,
            gex_library_metrics_summary,
            gex_bc_plot,
            gex_genes_plot,
            gex_sequencing_plot,
            crm_version,
        )

        _build_vdj_t_data(
            data["data"],
            sample_name,
            vdj_t_alerts,
            vdj_t_expression_metrics,
            vdj_t_annotation_metrics,
            vdj_t_bc_plot,
            crm_version,
        )

        _build_vdj_b_data(
            data["data"],
            sample_name,
            vdj_b_alerts,
            vdj_b_expression_metrics,
            vdj_b_annotation_metrics,
            vdj_b_bc_plot,
            crm_version,
        )

        _build_antibody_data(
            data["data"],
            sample_name,
            antibody_alerts,
            antibody_expression_mapping_metrics,
            antibody_bc_plot,
            crm_version,
        )

        module.add_data_source(file, sample_name, module="cellranger", section="multi")

    general_data.data, general_data.headers = combine_data_and_headers(
        [
            (
                gex_cells_metrics.data,
                gex_cells_metrics.headers,
                ["Cells", "Mean reads per Cell", "Median reads per cell", "Median genes per cell"],
            ),
            (
                gex_library_metrics_summary.data,
                gex_library_metrics_summary.headers,
                ["Sequencing saturation"],
            ),
        ]
    )

    general_data.data = module.ignore_samples(general_data.data)
    if len(general_data.data) == 0:
        return 0
    module.general_stats_addcols(general_data.data, general_data.headers)

    module.write_data_file(gex_cells_metrics.data, "multiqc_cellranger_gex_cells")
    module.write_data_file(gex_library_metrics_summary.data, "multiqc_cellranger_gex_libraries")
    module.write_data_file(vdj_t_expression_metrics.data, "multiqc_cellranger_vdj_t_expr")
    module.write_data_file(vdj_t_annotation_metrics.data, "multiqc_cellranger_vdj_t_anno")
    module.write_data_file(vdj_b_expression_metrics.data, "multiqc_cellranger_vdj_b_expr")
    module.write_data_file(vdj_b_annotation_metrics.data, "multiqc_cellranger_vdj_b_anno")
    module.write_data_file(antibody_expression_mapping_metrics.data, "multiqc_cellranger_antibody_expr")

    for section in [
        gex_alerts,
        gex_cells_metrics,
        gex_library_metrics_summary,
        gex_bc_plot,
        gex_genes_plot,
        gex_sequencing_plot,
        vdj_t_alerts,
        vdj_t_expression_metrics,
        vdj_t_annotation_metrics,
        vdj_t_bc_plot,
        vdj_b_alerts,
        vdj_b_expression_metrics,
        vdj_b_annotation_metrics,
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


def _get_column_colors(tab: str, section: str, version: str) -> Dict:
    is_v7 = re.match(r"7\.[0-9]\.[0-9]", version)
    is_v8 = re.match(r"8\.[0-9]\.[0-9]", version)
    if not (is_v7 or is_v8):
        log.warning(f"Version {version} did not match a supported CellRanger version")

    if tab == "gex":
        if section == "hero_metrics":
            if is_v7:
                return {
                    "Cells": "Purples",
                    "Median reads per cell": "Greens",
                    "Median genes per cell": "Blues",
                    "Total genes detected": "Purples",
                    "Median UMI counts per cell": "Greys",
                    "Confidently mapped reads in cells": "Greens",
                }
            if is_v8:
                return {
                    "Cells": "Purples",
                    "Mean reads per cell": "Greens",
                    "Median genes per cell": "Blues",
                    "Total genes detected": "Purples",
                    "Median UMI counts per cell": "Greys",
                    "Confidently mapped reads in cells": "Greens",
                }
            # Return mix of everything as last resort
            return {
                "Cells": "Purples",
                "Mean reads per cell": "Greens",
                "Median reads per cell": "Greens",
                "Median genes per cell": "Blues",
                "Total genes detected": "Purples",
                "Median UMI counts per cell": "Greys",
                "Confidently mapped reads in cells": "Greens",
            }
        if section == "sequencing_metrics_table":
            # Is the same for v7.* and v8.*
            return {
                "Fastq ID": "Greys",
                "Number of reads": "Greys",
                "Number of short reads skipped": "Greys",
                "Q30 barcodes": "Greens",
                "Q30 UMI": "Blues",
                "Q30 RNA read": "Greens",
                "Q30 RNA read 2": "Blues",
            }
        if section == "physical_library_metrics_table":
            # Is the same for v7.* and v8.*
            return {
                "Physical library ID": "Greys",
                "Number of reads": "Greys",
                "Valid barcodes": "Greens",
                "Valid UMIs": "Blues",
                "Sequencing saturation": "YlGn",
                "Confidently mapped reads in cells": "Greys",
                "Mean reads per cell": "Greys",
            }
        log.warning(f"Section {section} not recognized in tab {tab}")
        return {}
    if tab == "vdj-t":
        if section == "hero_metrics":
            return {
                "Estimated number of cells": "Purples",
                "Number of cells with productive V-J spanning pair": "Greens",
                "Median TRA UMIs per Cell": "Blues",
                "Median TRB UMIs per Cell": "Greens",
            }
        if section == "annotation_metrics_table":
            return {
                "Cells with productive V-J spanning pair": "Greens",
                "Cells with productive V-J spanning (TRA, TRB) pair": "Greys",  # hidden
                "Cells with productive TRA contig": "Blues",
                "Cells with productive TRB contig": "Greens",
                "Paired clonotype diversity": "Purples",
            }
        log.warning(f"Section {section} not recognized in tab {tab}")
        return {}
    if tab == "vdj-b":
        if section == "hero_metrics":
            return {
                "Estimated number of cells": "Purples",
                "Number of cells with productive V-J spanning pair": "Greens",
                "Median IGH UMIs per Cell": "Blues",
                "Median IGK UMIs per Cell": "Greens",
            }
        if section == "annotation_metrics_table":
            return {
                "Cells with productive V-J spanning pair": "Greens",
                "Cells with productive V-J spanning (IGK, IGH) pair": "Greys",  # hidden
                "Cells with productive V-J spanning (IGL, IGH) pair": "Greys",  # hidden
                "Cells with productive IGH contig": "Blues",
                "Cells with productive IGK contig": "Greens",
                "Cells with productive IGL contig": "Blues",
                "Paired clonotype diversity": "Purples",
            }
        log.warning(f"Section {section} not recognized in tab {tab}")
        return {}
    if tab == "antibody":
        if section == "hero_metrics":
            if is_v7:
                return {
                    "Cells": "Purples",
                    "Median UMI counts per cell": "Greens",
                    "Mean antibody reads usable per cell": "Blues",
                }
            if is_v8:
                return {
                    "Cells": "Purples",
                    "Median UMI counts per cell": "Greens",
                    "Mean antibody reads usable per cell": "Blues",
                    "Antibody reads in cells": "Greens",
                }
            return {
                "Cells": "Purples",
                "Median UMI counts per cell": "Greens",
                "Mean antibody reads usable per cell": "Blues",
                "Antibody reads in cells": "Greens",
            }
        log.warning(f"Section {section} not recognized in tab {tab}")
        return {}
    log.warning(f"Tab {tab} not recognized")
    return {}


def _parse_table(
    table_root: Dict, color_dict: Dict[str, str], filter_with_color_dict: bool = True
) -> Tuple[Dict, Dict]:
    if not table_root:
        return {}, {}

    if len(table_root["table"]["rows"]) > 1:
        log.error(f"Multiple rows for table: {table_root['help'].get('title', 'UNDEFINED TABLE TITLE')}")

    table_keys = table_root["table"]["header"]
    table_data = table_root["table"]["rows"][0]
    table_help = table_root["help"]

    parsed_data = dict(zip(table_keys, table_data))
    parsed_headers = {entry[0]: {"title": entry[0], "description": entry[1][0]} for entry in table_help["data"]}

    # if filter_with_color_dict is True then the keys in the color dict are used as filters for columns we want to add
    if filter_with_color_dict:
        filtered_data, filtered_headers = {}, {}
        for key in color_dict.keys():
            if key in parsed_data and key in parsed_headers:
                filtered_data[key] = parsed_data[key]
                filtered_headers[key] = parsed_headers[key]
            else:
                log.debug(f"Removing unexpected column {key}")
        parsed_data, parsed_headers = filtered_data, filtered_headers

    parsed_data, parsed_headers = update_data_and_header(parsed_data, parsed_headers, color_dict, log)
    return parsed_data, parsed_headers


def _parse_all_alerts(sample_alert_list: List, library_alert_list: List):
    def parse_alert_list(alert_list: List):
        def modify(level: str):
            if level == "INFO":
                return "PASS"
            if level == "ERROR":
                return "FAIL"
            return level

        return (
            {alert["title"]: alert["level"] for alert in alert_list},
            {alert["title"]: {"modify": modify} for alert in alert_list},
        )

    sample_alerts, sample_headers = parse_alert_list(sample_alert_list)
    library_alerts, library_headers = parse_alert_list(library_alert_list)
    total_alerts, total_headers = (library_alerts, library_headers)

    # Careful dictionary update:
    for sample_alert_key in sample_alerts.keys():
        if sample_alert_key not in total_alerts:
            total_alerts[sample_alert_key] = sample_alerts[sample_alert_key]
            total_headers[sample_alert_key] = sample_headers[sample_alert_key]
        else:
            # ERROR > WARN > INFO
            remap = {"ERROR": 2, "WARN": 1, "INFO": 0}
            if remap.get(sample_alerts[sample_alert_key], 3) > remap.get(total_alerts[sample_alert_key], 3):
                total_alerts[sample_alert_key] = sample_alerts[sample_alert_key]
                total_headers[sample_alert_key] = sample_headers[sample_alert_key]
    return total_alerts, total_headers


def _parse_plot(data: Dict):
    if not data:
        return None, None
    if ("plot" not in data) or ("data" not in data["plot"]):
        return None, None
    if ("help" not in data) or ("helpText" not in data["help"]):
        return None, None

    plot_data = {}
    for subset in data["plot"]["data"]:
        plot_data.update(dict(zip(subset["x"], subset["y"])))

    return plot_data, {"helptext": data["help"]["helpText"]}


def _build_gex_data(
    data: Dict,
    sample_name: str,
    alerts: SectionData,
    cells_metrics: SectionData,
    library_metrics: SectionData,
    barcode_plot: SectionData,
    genes_plot: SectionData,
    sequencing_plot: SectionData,
    version: str,
):
    sample_websummary = data.get("sample_websummary", {}).get("gex_tab", {})
    library_websummary = data.get("library_websummary", {}).get("gex_tab", {})

    if not sample_websummary and not library_websummary:
        log.debug("Gene-Expression data incomplete")
        return
    if not sample_websummary or not library_websummary:
        log.debug("No Gene-Expression data found")
        return

    parsed_alerts, parsed_alert_headers = _parse_all_alerts(sample_websummary["alerts"], library_websummary["alerts"])
    if parsed_alerts:
        alerts.update_sample(parsed_alerts, parsed_alert_headers, sample_name)

    try:
        cells_source = sample_websummary["content"]["hero_metrics"]
        cells, cells_headers = _parse_table(cells_source, _get_column_colors("gex", "hero_metrics", version))

        if cells:
            cells_metrics.update_sample(cells, cells_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/gex_tab/content/hero_metrics")

    try:
        sequencing_metrics_source = library_websummary["content"][
            "sequencing_metrics_table"
        ]  # could have multiple rows?
        sequencing_metrics, sequencing_metrics_headers = _parse_table(
            sequencing_metrics_source, _get_column_colors("gex", "sequencing_metrics_table", version)
        )

        if sequencing_metrics:
            library_metrics.update_sample(sequencing_metrics, sequencing_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/gex_tab/content/sequencing_metrics_table")

    try:
        physical_library_metrics_source = library_websummary["content"]["physical_library_metrics_table"]
        physical_library_metrics, physical_library_metrics_headers = _parse_table(
            physical_library_metrics_source, _get_column_colors("gex", "physical_library_metrics_table", version)
        )

        if physical_library_metrics:
            library_metrics.update_sample(physical_library_metrics, physical_library_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/gex_tab/content/physical_library_metrics_table")

    try:
        bc_plot_source: Dict = library_websummary.get("content", {}).get("barcode_rank_plot", {})
        if bc_plot_source:
            bc_plot_data, bc_plot_headers = _parse_plot(bc_plot_source)
            if bc_plot_data:
                barcode_plot.update_sample(bc_plot_data, bc_plot_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/gex_tab/content/barcode_rank_plot")

    try:
        genes_plot_source: Dict = library_websummary.get("content", {}).get("median_genes_per_cell_plot", {})
        if genes_plot_source:
            genes_plot_data, genes_plot_headers = _parse_plot(genes_plot_source)
            if genes_plot_data:
                genes_plot.update_sample(genes_plot_data, genes_plot_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/gex_tab/content/median_genes_per_cell_plot")

    try:
        sequencing_plot_source: Dict = library_websummary.get("content", {}).get("sequencing_saturation_plot", {})
        if sequencing_plot_source:
            sequencing_plot_data, sequencing_plot_headers = _parse_plot(sequencing_plot_source)
            if sequencing_plot_data:
                sequencing_plot.update_sample(sequencing_plot_data, sequencing_plot_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/gex_tab/content/sequencing_saturation_plot")


def _build_vdj_t_data(
    data: Dict,
    sample_name: str,
    alerts: SectionData,
    expression_metrics: SectionData,
    annotation_metrics: SectionData,
    barcode_plot: SectionData,
    version: str,
):
    sample_websummary = data.get("sample_websummary", {}).get("vdj_t_tab", {})
    library_websummary = data.get("library_websummary", {}).get("vdj_t_tab", {})

    if not sample_websummary and not library_websummary:
        log.debug("No VDJ-T data found")
        return
    if not sample_websummary or not library_websummary:
        log.debug("VDJ-T data incomplete")
        return

    parsed_alerts, parsed_alert_headers = _parse_all_alerts(sample_websummary["alerts"], library_websummary["alerts"])
    if parsed_alerts:
        alerts.update_sample(parsed_alerts, parsed_alert_headers, sample_name)

    try:
        expr_metrics_source = sample_websummary["content"]["hero_metrics"]
        expr_metrics, expr_metrics_headers = _parse_table(
            expr_metrics_source, _get_column_colors("vdj-t", "hero_metrics", version)
        )

        if expr_metrics:
            expression_metrics.update_sample(expr_metrics, expr_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/vdj_t_tab/content/hero_metrics")

    try:
        anno_metrics_source = sample_websummary["content"]["annotation_metrics_table"]
        anno_metrics, anno_metrics_headers = _parse_table(
            anno_metrics_source, _get_column_colors("vdj-t", "annotation_metrics_table", version)
        )

        if anno_metrics:
            annotation_metrics.update_sample(anno_metrics, anno_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/vdj_t_tab/content/annotation_metrics_table")

    try:
        bc_plot_source: Dict = library_websummary.get("content", {}).get("barcode_rank_plot", {})
        if bc_plot_source:
            bc_plot_data, bc_plot_headers = _parse_plot(bc_plot_source)
            if bc_plot_data:
                barcode_plot.update_sample(bc_plot_data, bc_plot_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/vdj_t_tab/content/barcode_rank_plot")


def _build_vdj_b_data(
    data: Dict,
    sample_name: str,
    alerts: SectionData,
    expression_metrics: SectionData,
    annotation_metrics: SectionData,
    barcode_plot: SectionData,
    version: str,
):
    sample_websummary = data.get("sample_websummary", {}).get("vdj_b_tab", {})
    library_websummary = data.get("library_websummary", {}).get("vdj_b_tab", {})

    if not sample_websummary and not library_websummary:
        log.debug("No VDJ-B data found")
        return
    if not sample_websummary or not library_websummary:
        log.debug("VDJ-B data incomplete")
        return

    parsed_alerts, parsed_alert_headers = _parse_all_alerts(sample_websummary["alerts"], library_websummary["alerts"])
    if parsed_alerts:
        alerts.update_sample(parsed_alerts, parsed_alert_headers, sample_name)

    try:
        expr_metrics_source = sample_websummary["content"]["hero_metrics"]
        expr_metrics, expr_metrics_headers = _parse_table(
            expr_metrics_source, _get_column_colors("vdj-b", "hero_metrics", version)
        )

        if expr_metrics:
            expression_metrics.update_sample(expr_metrics, expr_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/vdj_b_tab/content/hero_metrics")

    try:
        anno_metrics_source = sample_websummary["content"]["annotation_metrics_table"]
        anno_metrics, anno_metrics_headers = _parse_table(
            anno_metrics_source, _get_column_colors("vdj-b", "annotation_metrics_table", version)
        )

        if anno_metrics:
            annotation_metrics.update_sample(anno_metrics, anno_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/vdj_b_tab/content/annotation_metrics_table")

    try:
        bc_plot_source: Dict = library_websummary.get("content", {}).get("barcode_rank_plot", {})
        if bc_plot_source:
            bc_plot_data, bc_plot_headers = _parse_plot(bc_plot_source)
            if bc_plot_data:
                barcode_plot.update_sample(bc_plot_data, bc_plot_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/vdj_b_tab/content/barcode_rank_plot")


def _build_antibody_data(
    data: Dict,
    sample_name: str,
    alerts: SectionData,
    expression_mapping_metrics: SectionData,
    barcode_plot: SectionData,
    version: str,
):
    sample_websummary = data.get("sample_websummary", {}).get("antibody_tab", {})
    library_websummary = data.get("library_websummary", {}).get("antibody_tab", {})

    if not sample_websummary and not library_websummary:
        log.debug("No antibody data found")
        return
    if not sample_websummary or not library_websummary:
        log.debug("Antibody data incomplete")
        return

    parsed_alerts, parsed_alert_headers = _parse_all_alerts(sample_websummary["alerts"], library_websummary["alerts"])
    if parsed_alerts:
        alerts.update_sample(parsed_alerts, parsed_alert_headers, sample_name)

    try:
        expression_metrics_source = sample_websummary["content"]["hero_metrics"]
        expression_metrics, expression_metrics_headers = _parse_table(
            expression_metrics_source, _get_column_colors("antibody", "hero_metrics", version)
        )

        if expression_metrics:
            expression_mapping_metrics.update_sample(expression_metrics, expression_metrics_headers, sample_name)
    except KeyError:
        log.debug("Could not find sample_websummary/antibody_tab/content/hero_metrics")

    try:
        bc_plot_source: Dict = library_websummary.get("content", {}).get("barcode_rank_plot", {})
        if bc_plot_source:
            bc_plot_data, bc_plot_headers = _parse_plot(bc_plot_source)
            if bc_plot_data:
                barcode_plot.update_sample(bc_plot_data, bc_plot_headers, sample_name)
    except KeyError:
        log.debug("Could not find library_websummary/antibody_tab/content/barcode_rank_plot")
