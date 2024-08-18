"""MParse output from Cell Ranger count"""

import json
import logging
import re
from typing import Dict, Optional

from multiqc import config, BaseMultiqcModule
from multiqc.modules.cellranger.utils import set_hidden_cols, update_dict, parse_bcknee_data, transform_data
from multiqc.plots import bargraph, linegraph, table


log = logging.getLogger(__name__)


def parse_count_html(module: BaseMultiqcModule) -> int:
    """
    Cell Ranger count report parser
    """

    data_by_sample: Dict[str, Dict] = dict()
    antibody_data_by_sample: Dict[str, Dict] = dict()
    general_data_by_sample: Dict[str, Dict] = dict()
    warnings_by_sample: Dict[str, Dict] = dict()
    all_plots_params_by_id: Dict[str, Dict] = {"bc": dict(), "genes": dict()}
    all_plots_data_by_sname_by_id: Dict[str, Dict] = {"bc": dict(), "genes": dict()}
    general_data_headers: Dict[str, Dict] = dict()
    count_data_headers: Dict[str, Dict] = dict()
    antibody_data_headers: Dict[str, Dict] = dict()
    count_warnings_headers: Dict[str, Dict] = dict()

    for f in module.find_log_files("cellranger/count_html", filehandles=True):
        summary: Optional[Dict] = None
        for line in f["f"]:
            line = line.strip()
            if line.startswith("const data"):
                line = line.replace("const data = ", "")
                summary = json.loads(line)["summary"]
                break
        if summary is None:
            continue

        s_name = module.clean_s_name(summary["sample"]["id"], f)

        # Extract software version
        try:
            version_pair = summary["summary_tab"]["pipeline_info_table"]["rows"][-1]
            assert version_pair[0] == "Pipeline Version"
            version_match = re.search(r"cellranger-([\d\.]+)", version_pair[1])
            if version_match:
                module.add_software_version(version_match.group(1), s_name)
        except (KeyError, AssertionError):
            log.debug(f"Unable to parse version for sample {s_name}")

        data_general_stats: Dict[str, Dict] = dict()

        # Store general stats from cells
        col_dict = {
            "Estimated Number of Cells": "estimated cells",
            "Mean Reads per Cell": "avg reads/cell",
            "Fraction Reads in Cells": "reads in cells",
        }
        colours = {
            "estimated cells": "PuBu",
            "avg reads/cell": "GnBu",
            "reads in cells": "Purples",
        }
        data_general_stats, general_data_headers = update_dict(
            data_general_stats,
            general_data_headers,
            summary["summary_tab"]["cells"]["table"]["rows"],
            col_dict,
            colours,
            "Count",
        )

        # Store general stats from sequencing tables
        col_dict = {
            "Number of Reads": "reads",
            "Valid Barcodes": "valid bc",
            "Q30 Bases in Barcode": "Q30 bc",
            "Q30 Bases in UMI": "Q30 UMI",
            "Q30 Bases in RNA Read": "Q30 read",
        }
        colours = {
            "reads": "PuBuGn",
            "valid bc": "RdYlGn",
            "Q30 bc": "RdYlBu",
            "Q30 UMI": "Spectral",
            "Q30 read": "RdBu",
        }
        data_general_stats, general_data_headers = update_dict(
            data_general_stats,
            general_data_headers,
            summary["summary_tab"]["sequencing"]["table"]["rows"],
            col_dict,
            colours,
            "Count",
        )

        # Store full data from cell ranger count report
        data_rows = (
            summary["summary_tab"]["sequencing"]["table"]["rows"]
            + summary["summary_tab"]["cells"]["table"]["rows"]
            + summary["summary_tab"]["mapping"]["table"]["rows"]
        )
        col_dict = {
            "Number of Reads": "reads",
            "Estimated Number of Cells": "estimated cells",
            "Mean Reads per Cell": "avg reads/cell",
            "Total Genes Detected": "genes detected",
            "Median Genes per Cell": "median genes/cell",
            "Fraction Reads in Cells": "reads in cells",
            "Valid Barcodes": "valid bc",
            "Valid UMIs": "valid umi",
            "Median UMI Counts per Cell": "median umi/cell",
            "Sequencing Saturation": "saturation",
            "Q30 Bases in Barcode": "Q30 bc",
            "Q30 Bases in UMI": "Q30 UMI",
            "Q30 Bases in RNA Read": "Q30 read",
            "Reads Mapped to Genome": "reads mapped",
            "Reads Mapped Confidently to Genome": "confident reads",
            "Reads Mapped Confidently to Transcriptome": "confident transcriptome",
            "Reads Mapped Confidently to Exonic Regions": "confident exonic",
            "Reads Mapped Confidently to Intronic Regions": "confident intronic",
            "Reads Mapped Confidently to Intergenic Regions": "confident intergenic",
            "Reads Mapped Antisense to Gene": "reads antisense",
        }
        colours = {
            "reads": "YlGn",
            "estimated cells": "RdPu",
            "avg reads/cell": "Blues",
            "genes detected": "Greens",
            "median genes/cell": "Purples",
            "reads in cells": "PuBuGn",
            "valid bc": "Spectral",
            "valid umi": "RdYlGn",
            "median umi/cell": "YlGn",
            "saturation": "YlOrRd",
        }
        sample_data, count_data_headers = update_dict(
            data_general_stats,
            count_data_headers,
            data_rows,
            col_dict,
            colours,
            "Count",
        )
        if not sample_data:
            continue

        # Extract warnings if any
        warnings = dict()
        alarms_list = summary["alarms"].get("alarms", [])
        for alarm in alarms_list:
            # "Intron mode used" alarm added in Cell Ranger 7.0 lacks id
            if "id" not in alarm:
                continue
            warnings[alarm["id"]] = "FAIL"
            count_warnings_headers[alarm["id"]] = {
                "title": alarm["id"],
                "description": alarm["title"],
                "bgcols": {"FAIL": "#f06807"},
            }

        # Extract data for plots
        help_dict = {x[0]: x[1][0] for x in summary["summary_tab"]["cells"]["help"]["data"]}
        plots_params_by_id = {}
        plots_data_by_sname_by_id = {}
        if "analysis_tab" in summary:
            plots_params_by_id.update(
                {
                    "bc": {
                        "config": {
                            "id": "mqc_cellranger_count_bc_knee",
                            "title": f"Cell Ranger count: {summary['summary_tab']['cells']['barcode_knee_plot']['layout']['title']}",
                            "xlab": summary["summary_tab"]["cells"]["barcode_knee_plot"]["layout"]["xaxis"]["title"],
                            "ylab": summary["summary_tab"]["cells"]["barcode_knee_plot"]["layout"]["yaxis"]["title"],
                            "ylog": True,
                            "xlog": True,
                        },
                        "description": "Barcode knee plot",
                        "helptext": help_dict["Barcode Rank Plot"],
                    },
                    "genes": {
                        "config": {
                            "id": "mqc_cellranger_count_genesXcell",
                            "title": f"Cell Ranger count: {summary['analysis_tab']['median_gene_plot']['help']['title']}",
                            "xlab": summary["analysis_tab"]["median_gene_plot"]["plot"]["layout"]["xaxis"]["title"],
                            "ylab": summary["analysis_tab"]["median_gene_plot"]["plot"]["layout"]["yaxis"]["title"],
                            "ylog": False,
                            "xlog": False,
                        },
                        "description": "Median gene counts per cell",
                        "helptext": summary["analysis_tab"]["median_gene_plot"]["help"]["helpText"],
                    },
                }
            )
            try:
                plots_params_by_id["saturation"] = {
                    "config": {
                        "id": "mqc_cellranger_count_saturation",
                        "title": f"Cell Ranger count: {summary['analysis_tab']['seq_saturation_plot']['help']['title']}",
                        "xlab": summary["analysis_tab"]["seq_saturation_plot"]["plot"]["layout"]["xaxis"]["title"],
                        "ylab": summary["analysis_tab"]["seq_saturation_plot"]["plot"]["layout"]["yaxis"]["title"],
                        "ylog": False,
                        "xlog": False,
                        "ymin": 0,
                        "ymax": 1,
                    },
                    "description": "Sequencing saturation",
                    "helptext": summary["analysis_tab"]["seq_saturation_plot"]["help"]["helpText"],
                }
            except KeyError:
                pass

            plots_data_by_sname_by_id.update(
                {
                    "bc": parse_bcknee_data(summary["summary_tab"]["cells"]["barcode_knee_plot"]["data"], s_name),
                    "genes": {s_name: transform_data(summary["analysis_tab"]["median_gene_plot"]["plot"]["data"][0])},
                }
            )
            if "seq_saturation_plot" in summary["analysis_tab"]:
                plots_data_by_sname_by_id["saturation"] = {
                    s_name: transform_data(summary["analysis_tab"]["seq_saturation_plot"]["plot"]["data"][0])
                }

        # Store full data for ANTIBODY capture
        antibody_data: Dict[str, Dict] = dict()
        if "ANTIBODY_sequencing" in summary["summary_tab"]:
            data_rows = (
                summary["summary_tab"]["ANTIBODY_sequencing"]["table"]["rows"]
                + summary["summary_tab"]["ANTIBODY_application"]["table"]["rows"]
            )
            col_dict = {
                "Number of Reads": "reads",
                "Valid Barcodes": "valid bc",
                "Valid UMIs": "valid umi",
                "Sequencing Saturation": "saturation",
                "Q30 Bases in Barcode": "Q30 bc",
                "Q30 Bases in Antibody Read": "Q30 read",
                "Q30 Bases in UMI": "Q30 UMI",
                "Fraction Antibody Reads": "antibody reads",
                "Fraction Antibody Reads Usable": "antibody reads usable",
                "Antibody Reads Usable per Cell": "antibody reads usable/cell",
                "Fraction Antibody Reads in Aggregate Barcodes": "reads in aggregate bc",
                "Fraction Unrecognized Antibody": "unrecognized antibody",
                "Antibody Reads in Cells": "antibody reads in cells",
                "Median UMIs per Cell (summed over all recognized antibody barcodes)": "umi per cell",
            }
            colours = {
                "reads": "YlGn",
                "antibody reads": "RdPu",
                "reads in cells": "Blues",
                "reads usable": "Greens",
                "reads usable per cell": "Purples",
                "reads in aggregate bc": "PuBuGn",
                "valid bc": "Spectral",
                "valid umi": "RdYlGn",
                "Q30 bc": "YlGn",
                "saturation": "YlOrRd",
            }
            antibody_data, antibody_data_headers = update_dict(
                antibody_data,
                antibody_data_headers,
                data_rows,
                col_dict,
                colours,
                "Antibody",
            )

            # Extract labels and values for the bargraph data
            combined_data = {}
            if "antibody_treemap_plot" in summary["antibody_tab"]:
                for label, value in zip(
                    summary["antibody_tab"]["antibody_treemap_plot"]["plot"]["data"][0]["labels"],
                    summary["antibody_tab"]["antibody_treemap_plot"]["plot"]["data"][0]["values"],
                ):
                    label_match = re.search(r"<b>(.*?)\s+\((.*?)%\)</b>", label)
                    if label_match:
                        label_value = label_match.group(1)
                        value_ = round(value * 100, 2)
                        combined_data[label_value] = value_

                # Extract labels and number of cells for labelling the bargraph
                combined_label = {}
                for label, cells in zip(
                    summary["antibody_tab"]["antibody_treemap_plot"]["plot"]["data"][0]["labels"],
                    summary["antibody_tab"]["antibody_treemap_plot"]["plot"]["data"][0]["text"],
                ):
                    label_match = re.search(r"<b>(.*?)\s+\((.*?)%\)</b>", label)
                    if label_match:
                        label_value = label_match.group(1)
                        combined_label[label_value] = label_value + ": " + cells

                # Use the label from `combined_label` for the plot
                keys = dict()
                for key, value in combined_label.items():
                    keys[key] = {"name": value}

                plots_params_by_id["antibody_counts"] = {
                    "config": {
                        "id": "mqc_cellranger_antibody_counts",
                        "title": "Cell Ranger: Distribution of Antibody Counts",
                        "ylab": "% Total UMI",
                        "ymax": 100,
                        "cpswitch": False,
                        "use_legend": False,
                        "tt_decimals": 2,
                        "tt_suffix": "%",
                    },
                    "keys": keys,
                    "description": "Antibody Counts Distribution Plot",
                    "helptext": "Relative composition of antibody counts for features with at least 1 UMI. Box size represents fraction of total UMIs from cell barcodes that are derived from this antibody. Hover over a box to view more information on a particular antibody, including number of associated barcodes.",
                }
                plots_data_by_sname_by_id["antibody_counts"] = {s_name: combined_data}

        if s_name in general_data_by_sample:
            log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
        module.add_data_source(f, s_name, module="cellranger", section="count")
        data_by_sample[s_name] = sample_data
        if "antibody_tab" in summary:
            antibody_data_by_sample[s_name] = antibody_data
        general_data_by_sample[s_name] = data_general_stats
        if len(warnings) > 0:
            warnings_by_sample[s_name] = warnings
        all_plots_params_by_id.update(plots_params_by_id)
        for k in plots_data_by_sname_by_id.keys():
            if k not in all_plots_data_by_sname_by_id.keys():
                all_plots_data_by_sname_by_id[k] = dict()
            all_plots_data_by_sname_by_id[k].update(plots_data_by_sname_by_id[k])

    data_by_sample = module.ignore_samples(data_by_sample)
    if antibody_data_by_sample:
        antibody_data_by_sample = module.ignore_samples(antibody_data_by_sample)
    general_data_by_sample = module.ignore_samples(general_data_by_sample)
    warnings_by_sample = module.ignore_samples(warnings_by_sample)
    for k in all_plots_data_by_sname_by_id.keys():
        all_plots_data_by_sname_by_id[k] = module.ignore_samples(all_plots_data_by_sname_by_id[k])

    general_data_headers["reads"] = {
        "rid": "count_genstats_reads",
        "title": f"{config.read_count_prefix} Reads",
        "description": f"Number of reads ({config.read_count_desc})",
        "modify": lambda x: x * config.read_count_multiplier,
        "shared_key": "read_count",
        "namespace": "Count",
    }
    general_data_headers = set_hidden_cols(general_data_headers, ["Q30 bc", "Q30 UMI", "Q30 read"])

    count_data_headers["reads"] = {
        "rid": "count_data_reads",
        "title": f"{config.read_count_prefix} Reads",
        "description": f"Number of reads ({config.read_count_desc})",
        "modify": lambda x: x * config.read_count_multiplier,
    }
    count_data_headers = set_hidden_cols(
        count_data_headers,
        [
            "Q30 bc",
            "Q30 UMI",
            "Q30 read",
            "reads in cells",
            "avg reads/cell",
            "confident reads",
            "confident transcriptome",
            "confident intronic",
            "confident intergenic",
            "reads antisense",
            "saturation",
        ],
    )

    if antibody_data_by_sample:
        antibody_data_headers["reads"] = {
            "rid": "antibody_data_reads",
            "title": f"{config.read_count_prefix} Reads",
            "description": f"Number of reads ({config.read_count_desc})",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        antibody_data_headers = set_hidden_cols(
            antibody_data_headers,
            ["Q30 bc", "Q30 UMI", "Q30 read", "saturation", "umi per cell", "reads in aggregate bc"],
        )

    if len(general_data_by_sample) == 0:
        return 0

    module.general_stats_addcols(general_data_by_sample, general_data_headers)

    # Write parsed report data to a file
    module.write_data_file(data_by_sample, "multiqc_cellranger_count")
    if antibody_data_by_sample:
        module.write_data_file(antibody_data_by_sample, "multiqc_cellranger_antibody_count")

    # Add sections to the report
    if len(warnings_by_sample) > 0:
        module.add_section(
            name="Count - Warnings",
            anchor="cellranger-count-warnings",
            description="Warnings encountered during the analysis",
            plot=table.plot(
                warnings_by_sample,
                count_warnings_headers,
                {
                    "namespace": "Count",
                    "id": "cellranger-count-warnings-table",
                    "title": "Cellranger Count: Warnings",
                },
            ),
        )

    module.add_section(
        name="Count - Summary stats",
        anchor="cellranger-count-stats",
        description="Summary QC metrics from Cell Ranger count",
        plot=table.plot(
            data_by_sample,
            count_data_headers,
            {
                "namespace": "Count",
                "id": "cellranger-count-stats-table",
                "title": "Cellranger Count: Summary Stats",
            },
        ),
    )

    if all_plots_params_by_id.get("bc"):
        module.add_section(
            name="Count - BC rank plot",
            anchor="cellranger-count-bcrank-plot",
            description=all_plots_params_by_id["bc"].get("description"),
            helptext=all_plots_params_by_id["bc"]["helptext"],
            plot=linegraph.plot(all_plots_data_by_sname_by_id["bc"], all_plots_params_by_id["bc"]["config"]),
        )

    if all_plots_params_by_id.get("genes"):
        module.add_section(
            name="Count - Median genes",
            anchor="cellranger-count-genes-plot",
            description=all_plots_params_by_id["genes"]["description"],
            helptext=all_plots_params_by_id["genes"]["helptext"],
            plot=linegraph.plot(all_plots_data_by_sname_by_id["genes"], all_plots_params_by_id["genes"]["config"]),
        )

    if all_plots_params_by_id.get("saturation"):
        module.add_section(
            name="Count - Saturation plot",
            anchor="cellranger-count-saturation-plot",
            description=all_plots_params_by_id["saturation"]["description"],
            helptext=all_plots_params_by_id["saturation"]["helptext"],
            plot=linegraph.plot(
                all_plots_data_by_sname_by_id["saturation"],
                all_plots_params_by_id["saturation"]["config"],
            ),
        )

    if antibody_data_by_sample:
        module.add_section(
            name="Antibody - Summary stats",
            anchor="cellranger-antibody-stats",
            description="Summary QC metrics from Cell Ranger count",
            plot=table.plot(
                antibody_data_by_sample,
                antibody_data_headers,
                {
                    "namespace": "Antibody",
                    "id": "cellranger-antibody-stats-table",
                    "title": "Cellranger Antibody: Summary Stats",
                },
            ),
        )

    if all_plots_params_by_id.get("antibody_counts"):
        module.add_section(
            name="Antibody - Counts Distribution Bargraph",
            anchor="cellranger-antibody-counts",
            description=all_plots_params_by_id["antibody_counts"]["description"],
            helptext=all_plots_params_by_id["antibody_counts"]["helptext"],
            plot=bargraph.plot(
                all_plots_data_by_sname_by_id["antibody_counts"],
                all_plots_params_by_id["antibody_counts"]["keys"],
                all_plots_params_by_id["antibody_counts"]["config"],
            ),
        )

    return len(general_data_by_sample)
