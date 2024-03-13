""" MultiQC module to parse output from Cell Ranger count """

import json
import logging
import re

from multiqc import config
from multiqc.modules.cellranger.utils import set_hidden_cols, update_dict, parse_bcknee_data, clean_title_case
from multiqc.plots import linegraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class CellRangerVdjMixin:
    """Cell Ranger count report parser"""

    def parse_vdj_html(self):
        self.cellrangervdj_mapping = dict()
        self.cellrangervdj_annotations = dict()
        self.cellrangervdj_general_data = dict()
        self.cellrangervdj_warnings = dict()
        self.cellrangervdj_plots_conf = {"bc": dict(), "genes": dict()}
        self.cellrangervdj_plots_data = {"bc": dict(), "genes": dict()}
        self.vdj_general_data_headers = dict()
        self.vdj_mapping_headers = dict()
        self.vdj_annotations_headers = dict()
        self.vdj_warnings_headers = dict()

        for f in self.find_log_files("cellranger/vdj_html", filehandles=True):
            self.parse_vdj_report(f)

        self.cellrangervdj_mapping = self.ignore_samples(self.cellrangervdj_mapping)
        self.cellrangervdj_annotations = self.ignore_samples(self.cellrangervdj_annotations)
        self.cellrangervdj_general_data = self.ignore_samples(self.cellrangervdj_general_data)
        self.cellrangervdj_warnings = self.ignore_samples(self.cellrangervdj_warnings)
        for k in self.cellrangervdj_plots_data.keys():
            self.cellrangervdj_plots_data[k] = self.ignore_samples(self.cellrangervdj_plots_data[k])

        self.vdj_general_data_headers["reads"] = {
            "title": f"{config.read_count_prefix} Reads",
            "description": f"Number of reads ({config.read_count_desc})",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "namespace": "VDJ",
        }
        self.vdj_general_data_headers = set_hidden_cols(
            self.vdj_general_data_headers,
            [
                "Q30 bc",
                "Q30 UMI",
                "Q30 read",
                "Q30 read1",
                "Q30 read2",
            ],
        )

        self.vdj_mapping_headers["reads"] = {
            "title": f"{config.read_count_prefix} Reads",
            "description": f"Number of reads ({config.read_count_desc})",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.vdj_mapping_headers = set_hidden_cols(
            self.vdj_mapping_headers,
            [
                "Q30 bc",
                "Q30 UMI",
                "Q30 read1",
                "Q30 read2",
                "IGH reads",
                "IGK reads",
                "IGL reads",
            ],
        )

        self.vdj_annotations_headers = set_hidden_cols(
            self.vdj_annotations_headers,
            [
                "cells IGH contig",
                "cells IGK contig",
                "cells IGL contig",
                "cells IGH CDR3",
                "cells IGK CDR3",
                "cells IGL CDR3",
                "cells IGH VJ span",
                "cells IGK VJ span",
                "cells IGL VJ span",
            ],
        )

        if len(self.cellrangervdj_general_data) == 0:
            return 0

        else:
            for k in self.vdj_general_data_headers.keys():
                self.vdj_general_data_headers[k]["title"] = f"{self.vdj_general_data_headers[k]['title']} (VDJ)"
                self.vdj_general_data_headers[k][
                    "description"
                ] = f"{self.vdj_general_data_headers[k]['description']} (VDJ)"

            self.general_stats_addcols(self.cellrangervdj_general_data, self.vdj_general_data_headers)

            # Write parsed report data to a file
            self.write_data_file(self.cellrangervdj_mapping, "multiqc_cellranger_vdj_mapping")
            self.write_data_file(self.cellrangervdj_annotations, "multiqc_cellranger_vdj_annotations")

            # Add sections to the report
            if len(self.cellrangervdj_warnings) > 0:
                self.add_section(
                    name="VDJ - Warnings",
                    anchor="cellranger-vdj-warnings",
                    description="Warnings encountered during the analysis",
                    plot=table.plot(
                        self.cellrangervdj_warnings,
                        self.vdj_warnings_headers,
                        {
                            "namespace": "VDJ",
                            "id": "cellranger-vdj-warnings-table",
                            "title": "Cellranger VDJ: Warnings",
                        },
                    ),
                )

            self.add_section(
                name="VDJ - Summary stats",
                anchor="cellranger-vdj-stats",
                description="Summary QC metrics from Cell Ranger count",
                plot=table.plot(
                    self.cellrangervdj_mapping,
                    self.vdj_mapping_headers,
                    {
                        "namespace": "VDJ",
                        "id": "cellranger-vdj-stats-table",
                        "title": "Cellranger VDJ: Summary stats",
                    },
                ),
            )

            self.add_section(
                name="VDJ - Annotations",
                anchor="cellranger-vdj-annot",
                description="V(D)J annotations from Cell Ranger VDJ analysis",
                plot=table.plot(
                    self.cellrangervdj_annotations,
                    self.vdj_annotations_headers,
                    {
                        "namespace": "VDJ",
                        "id": "cellranger-vdj-annot-table",
                        "title": "Cellranger VDJ: Annotations",
                    },
                ),
            )

            self.add_section(
                name="VDJ - BC rank plot",
                anchor="cellranger-vdj-bcrank-plot",
                description=self.cellrangervdj_plots_conf["bc"]["description"],
                helptext=self.cellrangervdj_plots_conf["bc"]["helptext"],
                plot=linegraph.plot(self.cellrangervdj_plots_data["bc"], self.cellrangervdj_plots_conf["bc"]["config"]),
            )

            return len(self.cellrangervdj_general_data)

    def parse_vdj_report(self, f):
        """Go through the html report of cell ranger and extract the data in a dicts"""

        for line in f["f"]:
            line = line.strip()
            if line.startswith("const data"):
                line = line.replace("const data = ", "")
                mydict = json.loads(line)
                mydict = mydict["summary"]
                break

        s_name = self.clean_s_name(mydict["sample"]["id"], f)

        # Extract software version
        try:
            version_pair = mydict["summary_tab"]["pipeline_info_table"]["rows"][-1]
            assert version_pair[0] == "Pipeline Version"
            version_match = re.search(r"cellranger-([\d\.]+)", version_pair[1])
            if version_match:
                self.add_software_version(version_match.group(1), s_name)
        except (KeyError, AssertionError):
            log.debug(f"Unable to parse version for sample {s_name}")

        data_general_stats = dict()

        # Store general stats from cells and sequencing tables
        col_dict = {
            "Estimated Number of Cells": "estimated cells",
            "Mean Reads per Cell": "avg reads/cell",
            "Fraction Reads in Cells": "reads in cells",
        }
        colours = {
            "estimated cells": "RdPu",
            "avg reads/cell": "Blues",
            "reads in cells": "PuBuGn",
        }
        data_general_stats, self.vdj_general_data_headers = update_dict(
            data_general_stats,
            self.vdj_general_data_headers,
            mydict["summary_tab"]["cells"]["table"]["rows"],
            col_dict,
            colours,
            "VDJ",
        )

        col_dict = {
            "Number of Reads": "reads",
            "Valid Barcodes": "valid bc",
            "Q30 Bases in Barcode": "Q30 bc",
            "Q30 Bases in UMI": "Q30 UMI",
            "Q30 Bases in RNA Read 1": "Q30 read1",
            "Q30 Bases in RNA Read 2": "Q30 read2",
        }
        data_general_stats, self.vdj_general_data_headers = update_dict(
            data_general_stats,
            self.vdj_general_data_headers,
            mydict["summary_tab"]["vdj_sequencing"]["table"]["rows"],
            col_dict,
            {},
            "VDJ",
        )

        # Store sequencing results data from vdj report
        data_rows = (
            mydict["summary_tab"]["vdj_sequencing"]["table"]["rows"]
            + mydict["summary_tab"]["cells"]["table"]["rows"]
            + mydict["summary_tab"]["vdj_enrichment"]["table"]["rows"]
        )
        col_dict = {
            "Number of Reads": "reads",
            "Valid Barcodes": "valid bc",
            "Estimated Number of Cells": "estimated cells",
            "Mean Reads per Cell": "avg reads/cell",
            "Mean Used Reads per Cell": "used reads/cell",
            "Fraction Reads in Cells": "reads in cells",
            "Reads Mapped to Any V(D)J Gene": "VDJ reads",
            "Reads Mapped to IGH": "IGH reads",
            "Reads Mapped to IGK": "IGK reads",
            "Reads Mapped to IGL": "IGL reads",
            "Q30 Bases in Barcode": "Q30 bc",
            "Q30 Bases in UMI": "Q30 UMI",
            "Q30 Bases in RNA Read 1": "Q30 read1",
            "Q30 Bases in RNA Read 2": "Q30 read2",
        }
        data, self.vdj_mapping_headers = update_dict(
            data_general_stats,
            self.vdj_mapping_headers,
            data_rows,
            col_dict,
            {},
            "VDJ",
        )

        # Store VDJ annotation and expression data from vdj report
        col_dict = {
            "Cells With Productive V-J Spanning Pair": "cells VJ span",
            "Cells With Productive V-J Spanning (IGK, IGH) Pair": "cells IGK-IGH span",
            "Cells With Productive V-J Spanning (IGL, IGH) Pair": "cells IGL-IGH span",
            "Paired Clonotype Diversity": "clonotype diversity",
            "Cells With IGH Contig": "cells IGH contig",
            "Cells With IGK Contig": "cells IGK contig",
            "Cells With IGL Contig": "cells IGL contig",
            "Cells With CDR3-annotated IGH Contig": "cells IGH CDR3",
            "Cells With CDR3-annotated IGK Contig": "cells IGK CDR3",
            "Cells With CDR3-annotated IGL Contig": "cells IGL CDR3",
            "Cells With V-J Spanning IGH Contig": "cells IGH VJ span",
            "Cells With V-J Spanning IGK Contig": "cells IGK VJ span",
            "Cells With V-J Spanning IGL Contig": "cells IGL VJ span",
            "Cells With Productive IGH Contig": "cells IGH productive",
            "Cells With Productive IGK Contig": "cells IGK productive",
            "Cells With Productive IGL Contig": "cells IGL productive",
        }
        colours = {
            # "cells VJ span": "",
            # "cells IGK-IGH span": "",
            # "cells IGL-IGH span": "",
            # "clonotype diversity": "",
            # "cells IGH contig": "",
            # "cells IGK contig": "",
            # "cells IGL contig": "",
            # "cells IGH CDR3": "",
            # "cells IGK CDR3": "",
            # "cells IGL CDR3": "",
            # "cells IGH VJ span": "",
            # "cells IGK VJ span": "",
            # "cells IGL VJ span": "",
            # "cells IGH productive": "",
            # "cells IGK productive": "",
            # "cells IGL productive": "",
        }
        data_annotations, self.vdj_annotations_headers = update_dict(
            data_general_stats,
            self.vdj_annotations_headers,
            mydict["summary_tab"]["vdj_annotation"]["table"]["rows"],
            col_dict,
            colours,
            "VDJ",
        )

        # Extract warnings if any
        warnings = dict()
        alarms_list = mydict["alarms"].get("alarms", [])
        for alarm in alarms_list:
            warnings[alarm["id"]] = "FAIL"
            self.vdj_warnings_headers[alarm["id"]] = {
                "title": clean_title_case(alarm["id"].replace("_", " ")),
                "description": alarm["title"],
                "bgcols": {"FAIL": "#f06807"},
            }

        # Extract data for plots
        help_dict = {x[0]: x[1][0] for x in mydict["summary_tab"]["cells"]["help"]["data"]}
        plots = {
            "bc": {
                "config": {
                    "id": "mqc_cellranger_vdj_bc_knee",
                    "title": f"Cell Ranger VDJ: {mydict['summary_tab']['cells']['barcode_knee_plot']['layout']['title']}",
                    "xlab": mydict["summary_tab"]["cells"]["barcode_knee_plot"]["layout"]["xaxis"]["title"],
                    "ylab": mydict["summary_tab"]["cells"]["barcode_knee_plot"]["layout"]["yaxis"]["title"],
                    "yLog": True,
                    "xLog": True,
                },
                "description": "Barcode knee plot",
                "helptext": help_dict["Barcode Rank Plot"],
            }
        }
        plots_data = {"bc": parse_bcknee_data(mydict["summary_tab"]["cells"]["barcode_knee_plot"]["data"], s_name)}

        if len(data) > 0:
            if s_name in self.cellrangervdj_general_data:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            self.add_data_source(f, s_name, module="cellranger", section="count")
            self.cellrangervdj_mapping[s_name] = data
            self.cellrangervdj_general_data[s_name] = data_general_stats
            self.cellrangervdj_annotations[s_name] = data_annotations
            if len(warnings) > 0:
                self.cellrangervdj_warnings[s_name] = warnings
            self.cellrangervdj_plots_conf = plots
            for k in plots_data.keys():
                if k not in self.cellrangervdj_plots_data.keys():
                    self.cellrangervdj_plots_data[k] = dict()
                self.cellrangervdj_plots_data[k].update(plots_data[k])
