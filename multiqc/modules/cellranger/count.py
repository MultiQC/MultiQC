""" MultiQC module to parse output from Cell Ranger count """

import json
import logging
import re

from multiqc import config
from multiqc.modules.cellranger.utils import set_hidden_cols, update_dict, parse_bcknee_data, transform_data
from multiqc.plots import bargraph, linegraph, table


# Initialise the logger
log = logging.getLogger(__name__)


class CellRangerCountMixin:
    """Cell Ranger count report parser"""

    def parse_count_html(self):
        self.cellrangercount_data = dict()
        self.cellrangercount_antibody_data = dict()
        self.cellrangercount_general_data = dict()
        self.cellrangercount_warnings = dict()
        self.cellrangercount_plots_conf = {"bc": dict(), "genes": dict()}
        self.cellrangercount_plots_data = {"bc": dict(), "genes": dict()}
        self.count_general_data_headers = dict()
        self.count_data_headers = dict()
        self.antibody_data_headers = dict()
        self.count_warnings_headers = dict()

        for f in self.find_log_files("cellranger/count_html", filehandles=True):
            self.parse_count_report(f)

        self.cellrangercount_data = self.ignore_samples(self.cellrangercount_data)
        if self.cellrangercount_antibody_data:
            self.cellrangercount_antibody_data = self.ignore_samples(self.cellrangercount_antibody_data)
        self.cellrangercount_general_data = self.ignore_samples(self.cellrangercount_general_data)
        self.cellrangercount_warnings = self.ignore_samples(self.cellrangercount_warnings)
        for k in self.cellrangercount_plots_data.keys():
            self.cellrangercount_plots_data[k] = self.ignore_samples(self.cellrangercount_plots_data[k])

        self.count_general_data_headers["reads"] = {
            "rid": "count_genstats_reads",
            "title": f"{config.read_count_prefix} Reads",
            "description": f"Number of reads ({config.read_count_desc})",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "namespace": "Count",
        }
        self.count_general_data_headers = set_hidden_cols(
            self.count_general_data_headers, ["Q30 bc", "Q30 UMI", "Q30 read"]
        )

        self.count_data_headers["reads"] = {
            "rid": "count_data_reads",
            "title": f"{config.read_count_prefix} Reads",
            "description": f"Number of reads ({config.read_count_desc})",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.count_data_headers = set_hidden_cols(
            self.count_data_headers,
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

        if self.cellrangercount_antibody_data:
            self.antibody_data_headers["reads"] = {
                "rid": "antibody_data_reads",
                "title": f"{config.read_count_prefix} Reads",
                "description": f"Number of reads ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
            }
            self.antibody_data_headers = set_hidden_cols(
                self.antibody_data_headers,
                ["Q30 bc", "Q30 UMI", "Q30 read", "saturation", "umi per cell", "reads in aggregate bc"],
            )

        if len(self.cellrangercount_general_data) == 0:
            return 0

        self.general_stats_addcols(self.cellrangercount_general_data, self.count_general_data_headers)

        # Write parsed report data to a file
        self.write_data_file(self.cellrangercount_data, "multiqc_cellranger_count")
        if self.cellrangercount_antibody_data:
            self.write_data_file(self.cellrangercount_antibody_data, "multiqc_cellranger_antibody_count")

        # Add sections to the report
        if len(self.cellrangercount_warnings) > 0:
            self.add_section(
                name="Count - Warnings",
                anchor="cellranger-count-warnings",
                description="Warnings encountered during the analysis",
                plot=table.plot(
                    self.cellrangercount_warnings,
                    self.count_warnings_headers,
                    {
                        "namespace": "Count",
                        "id": "cellranger-count-warnings-table",
                        "title": "Cellranger Count: Warnings",
                    },
                ),
            )

        self.add_section(
            name="Count - Summary stats",
            anchor="cellranger-count-stats",
            description="Summary QC metrics from Cell Ranger count",
            plot=table.plot(
                self.cellrangercount_data,
                self.count_data_headers,
                {
                    "namespace": "Count",
                    "id": "cellranger-count-stats-table",
                    "title": "Cellranger Count: Summary Stats",
                },
            ),
        )

        self.add_section(
            name="Count - BC rank plot",
            anchor="cellranger-count-bcrank-plot",
            description=self.cellrangercount_plots_conf["bc"]["description"],
            helptext=self.cellrangercount_plots_conf["bc"]["helptext"],
            plot=linegraph.plot(self.cellrangercount_plots_data["bc"], self.cellrangercount_plots_conf["bc"]["config"]),
        )

        self.add_section(
            name="Count - Median genes",
            anchor="cellranger-count-genes-plot",
            description=self.cellrangercount_plots_conf["genes"]["description"],
            helptext=self.cellrangercount_plots_conf["genes"]["helptext"],
            plot=linegraph.plot(
                self.cellrangercount_plots_data["genes"], self.cellrangercount_plots_conf["genes"]["config"]
            ),
        )

        if "saturation" in self.cellrangercount_plots_data:
            self.add_section(
                name="Count - Saturation plot",
                anchor="cellranger-count-saturation-plot",
                description=self.cellrangercount_plots_conf["saturation"]["description"],
                helptext=self.cellrangercount_plots_conf["saturation"]["helptext"],
                plot=linegraph.plot(
                    self.cellrangercount_plots_data["saturation"],
                    self.cellrangercount_plots_conf["saturation"]["config"],
                ),
            )

        if self.cellrangercount_antibody_data:
            self.add_section(
                name="Antibody - Summary stats",
                anchor="cellranger-antibody-stats",
                description="Summary QC metrics from Cell Ranger count",
                plot=table.plot(
                    self.cellrangercount_antibody_data,
                    self.antibody_data_headers,
                    {
                        "namespace": "Antibody",
                        "id": "cellranger-antibody-stats-table",
                        "title": "Cellranger Antibody: Summary Stats",
                    },
                ),
            )

        if "antibody_counts" in self.cellrangercount_plots_conf:
            self.add_section(
                name="Antibody - Counts Distribution Bargraph",
                anchor="cellranger-antibody-counts",
                description=self.cellrangercount_plots_conf["antibody_counts"]["description"],
                helptext=self.cellrangercount_plots_conf["antibody_counts"]["helptext"],
                plot=bargraph.plot(
                    self.cellrangercount_plots_data["antibody_counts"],
                    self.cellrangercount_plots_conf["antibody_counts"]["keys"],
                    self.cellrangercount_plots_conf["antibody_counts"]["config"],
                ),
            )

        return len(self.cellrangercount_general_data)

    def parse_count_report(self, f):
        """Go through the html report of cell ranger and extract the data in a dicts"""

        for line in f["f"]:
            line = line.strip()
            if line.startswith("const data"):
                line = line.replace("const data = ", "")
                summary = json.loads(line)
                summary = summary["summary"]
                break

        s_name = self.clean_s_name(summary["sample"]["id"], f)

        # Extract software version
        try:
            version_pair = summary["summary_tab"]["pipeline_info_table"]["rows"][-1]
            assert version_pair[0] == "Pipeline Version"
            version_match = re.search(r"cellranger-([\d\.]+)", version_pair[1])
            if version_match:
                self.add_software_version(version_match.group(1), s_name)
        except (KeyError, AssertionError):
            log.debug(f"Unable to parse version for sample {s_name}")

        data_general_stats = dict()

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
        data_general_stats, self.count_general_data_headers = update_dict(
            data_general_stats,
            self.count_general_data_headers,
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
        data_general_stats, self.count_general_data_headers = update_dict(
            data_general_stats,
            self.count_general_data_headers,
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
        table, self.count_data_headers = update_dict(
            data_general_stats,
            self.count_data_headers,
            data_rows,
            col_dict,
            colours,
            "Count",
        )
        if not table:
            return None

        # Extract warnings if any
        warnings = dict()
        alarms_list = summary["alarms"].get("alarms", [])
        for alarm in alarms_list:
            # "Intron mode used" alarm added in Cell Ranger 7.0 lacks id
            if "id" not in alarm:
                continue
            warnings[alarm["id"]] = "FAIL"
            self.count_warnings_headers[alarm["id"]] = {
                "title": alarm["id"],
                "description": alarm["title"],
                "bgcols": {"FAIL": "#f06807"},
            }

        # Extract data for plots
        help_dict = {x[0]: x[1][0] for x in summary["summary_tab"]["cells"]["help"]["data"]}
        plots = {
            "bc": {
                "config": {
                    "id": "mqc_cellranger_count_bc_knee",
                    "title": f"Cell Ranger count: {summary['summary_tab']['cells']['barcode_knee_plot']['layout']['title']}",
                    "xlab": summary["summary_tab"]["cells"]["barcode_knee_plot"]["layout"]["xaxis"]["title"],
                    "ylab": summary["summary_tab"]["cells"]["barcode_knee_plot"]["layout"]["yaxis"]["title"],
                    "yLog": True,
                    "xLog": True,
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
                    "yLog": False,
                    "xLog": False,
                },
                "description": "Median gene counts per cell",
                "helptext": summary["analysis_tab"]["median_gene_plot"]["help"]["helpText"],
            },
        }
        try:
            plots["saturation"] = {
                "config": {
                    "id": "mqc_cellranger_count_saturation",
                    "title": f"Cell Ranger count: {summary['analysis_tab']['seq_saturation_plot']['help']['title']}",
                    "xlab": summary["analysis_tab"]["seq_saturation_plot"]["plot"]["layout"]["xaxis"]["title"],
                    "ylab": summary["analysis_tab"]["seq_saturation_plot"]["plot"]["layout"]["yaxis"]["title"],
                    "yLog": False,
                    "xLog": False,
                    "ymin": 0,
                    "ymax": 1,
                },
                "description": "Sequencing saturation",
                "helptext": summary["analysis_tab"]["seq_saturation_plot"]["help"]["helpText"],
            }
        except KeyError:
            pass

        plots_data = {
            "bc": parse_bcknee_data(summary["summary_tab"]["cells"]["barcode_knee_plot"]["data"], s_name),
            "genes": {s_name: transform_data(summary["analysis_tab"]["median_gene_plot"]["plot"]["data"][0])},
        }
        if "seq_saturation_plot" in summary["analysis_tab"]:
            plots_data["saturation"] = {
                s_name: transform_data(summary["analysis_tab"]["seq_saturation_plot"]["plot"]["data"][0])
            }

        # Store full data for ANTIBODY capture
        antibody_data = dict()
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
            antibody_data, self.antibody_data_headers = update_dict(
                antibody_data,
                self.antibody_data_headers,
                data_rows,
                col_dict,
                colours,
                "Antibody",
            )

            # Extract labels and values for the bargraph data
            combined_data = {}
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

            plots["antibody_counts"] = {
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
            plots_data["antibody_counts"] = {s_name: combined_data}

        if s_name in self.cellrangercount_general_data:
            log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
        self.add_data_source(f, s_name, module="cellranger", section="count")
        self.cellrangercount_data[s_name] = table
        if "antibody_tab" in summary:
            self.cellrangercount_antibody_data[s_name] = antibody_data
        self.cellrangercount_general_data[s_name] = data_general_stats
        if len(warnings) > 0:
            self.cellrangercount_warnings[s_name] = warnings
        self.cellrangercount_plots_conf.update(plots)
        for k in plots_data.keys():
            if k not in self.cellrangercount_plots_data.keys():
                self.cellrangercount_plots_data[k] = dict()
            self.cellrangercount_plots_data[k].update(plots_data[k])
