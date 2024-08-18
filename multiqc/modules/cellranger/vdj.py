"""MultiQC module to parse output from Cell Ranger count"""

import json
import logging
import re
from typing import Dict

from multiqc import config, BaseMultiqcModule
from multiqc.modules.cellranger.utils import set_hidden_cols, update_dict, parse_bcknee_data, clean_title_case
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


def parse_vdj_html(module: BaseMultiqcModule) -> int:
    """
    Cell Ranger count report parser
    """
    mapping_by_sample: Dict[str, Dict] = dict()
    annotations_by_sample: Dict[str, Dict] = dict()
    general_data_by_sample: Dict[str, Dict] = dict()
    warnings_by_sample: Dict[str, Dict] = dict()
    plots_params_by_id: Dict[str, Dict] = {"bc": dict(), "genes": dict()}
    plots_data_by_sample_by_id: Dict[str, Dict] = {"bc": dict(), "genes": dict()}
    general_data_headers: Dict[str, Dict] = dict()
    mapping_headers: Dict[str, Dict] = dict()
    annotations_headers: Dict[str, Dict] = dict()
    vdj_warnings_headers: Dict[str, Dict] = dict()

    for f in module.find_log_files("cellranger/vdj_html", filehandles=True):
        mydict = None
        for line in f["f"]:
            line = line.strip()
            if line.startswith("const data"):
                line = line.replace("const data = ", "")
                mydict = json.loads(line)
                mydict = mydict["summary"]
                break
        if mydict is None:
            log.debug(f"Could not find VDJ data in {f['fn']}")
            continue
        if "vdj_sequencing" not in mydict["summary_tab"] or "cells" not in mydict["summary_tab"]:
            log.debug(f"Could not find VDJ sections in {f['fn']}")
            continue

        s_name = module.clean_s_name(mydict["sample"]["id"], f)

        # Extract software version
        try:
            version_pair = mydict["summary_tab"]["pipeline_info_table"]["rows"][-1]
            assert version_pair[0] == "Pipeline Version"
            version_match = re.search(r"cellranger-([\d\.]+)", version_pair[1])
            if version_match:
                module.add_software_version(version_match.group(1), s_name)
        except (KeyError, AssertionError):
            log.debug(f"Unable to parse version for sample {s_name}")

        data_general_stats: Dict[str, Dict] = dict()

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
        data_general_stats, general_data_headers = update_dict(
            data_general_stats,
            general_data_headers,
            mydict["summary_tab"]["cells"]["table"]["rows"],
            col_dict,
            colours,
            "VDJ",
        )

        # Store sequencing results data from vdj report
        col_dict = {
            "Number of Reads": "reads",
            "Valid Barcodes": "valid bc",
            "Q30 Bases in Barcode": "Q30 bc",
            "Q30 Bases in UMI": "Q30 UMI",
            "Q30 Bases in RNA Read 1": "Q30 read1",
            "Q30 Bases in RNA Read 2": "Q30 read2",
        }
        data_general_stats, general_data_headers = update_dict(
            data_general_stats,
            general_data_headers,
            mydict["summary_tab"]["vdj_sequencing"]["table"]["rows"],
            col_dict,
            {},
            "VDJ",
        )

        data_rows = mydict["summary_tab"]["vdj_sequencing"]["table"]["rows"]
        if "cells" in mydict["summary_tab"]:
            data_rows += mydict["summary_tab"]["cells"]["table"]["rows"]
        if "vdj_enrichment" in mydict["summary_tab"]:
            data_rows += mydict["summary_tab"]["vdj_enrichment"]["table"]["rows"]
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
        data, mapping_headers = update_dict(
            data_general_stats,
            mapping_headers,
            data_rows,
            col_dict,
            {},
            "VDJ",
        )

        if "vdj_annotation" in mydict["summary_tab"]:
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
            data_annotations, annotations_headers = update_dict(
                data_general_stats,
                annotations_headers,
                mydict["summary_tab"]["vdj_annotation"]["table"]["rows"],
                col_dict,
                {},
                "VDJ",
            )
        else:
            data_annotations = None

        # Extract warnings if any
        warnings = dict()
        alarms_list = mydict["alarms"].get("alarms", [])
        for alarm in alarms_list:
            warnings[alarm["id"]] = "FAIL"
            vdj_warnings_headers[alarm["id"]] = {
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
                    "ylog": True,
                    "xlog": True,
                },
                "description": "Barcode knee plot",
                "helptext": help_dict["Barcode Rank Plot"],
            }
        }
        plots_data = {"bc": parse_bcknee_data(mydict["summary_tab"]["cells"]["barcode_knee_plot"]["data"], s_name)}

        if len(data) > 0:
            if s_name in general_data_by_sample:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            module.add_data_source(f, s_name, module="cellranger", section="count")
            mapping_by_sample[s_name] = data
            general_data_by_sample[s_name] = data_general_stats
            if data_annotations:
                annotations_by_sample[s_name] = data_annotations
            if len(warnings) > 0:
                warnings_by_sample[s_name] = warnings
            plots_params_by_id = plots
            for k in plots_data.keys():
                if k not in plots_data_by_sample_by_id.keys():
                    plots_data_by_sample_by_id[k] = dict()
                plots_data_by_sample_by_id[k].update(plots_data[k])

    mapping_by_sample = module.ignore_samples(mapping_by_sample)
    annotations_by_sample = module.ignore_samples(annotations_by_sample)
    general_data_by_sample = module.ignore_samples(general_data_by_sample)
    warnings_by_sample = module.ignore_samples(warnings_by_sample)
    for k in plots_data_by_sample_by_id.keys():
        plots_data_by_sample_by_id[k] = module.ignore_samples(plots_data_by_sample_by_id[k])

    general_data_headers["reads"] = {
        "title": f"{config.read_count_prefix} Reads",
        "description": f"Number of reads ({config.read_count_desc})",
        "modify": lambda x: x * config.read_count_multiplier,
        "shared_key": "read_count",
        "namespace": "VDJ",
    }
    general_data_headers = set_hidden_cols(
        general_data_headers,
        [
            "Q30 bc",
            "Q30 UMI",
            "Q30 read",
            "Q30 read1",
            "Q30 read2",
        ],
    )

    mapping_headers["reads"] = {
        "title": f"{config.read_count_prefix} Reads",
        "description": f"Number of reads ({config.read_count_desc})",
        "modify": lambda x: x * config.read_count_multiplier,
    }
    mapping_headers = set_hidden_cols(
        mapping_headers,
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

    annotations_headers = set_hidden_cols(
        annotations_headers,
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

    if len(general_data_by_sample) == 0:
        return 0

    for k in general_data_headers.keys():
        general_data_headers[k]["title"] = f"{general_data_headers[k]['title']} (VDJ)"
        general_data_headers[k]["description"] = f"{general_data_headers[k]['description']} (VDJ)"

    module.general_stats_addcols(general_data_by_sample, general_data_headers)

    # Write parsed report data to a file
    module.write_data_file(mapping_by_sample, "multiqc_cellranger_vdj_mapping")
    module.write_data_file(annotations_by_sample, "multiqc_cellranger_vdj_annotations")

    # Add sections to the report
    if len(warnings_by_sample) > 0:
        module.add_section(
            name="VDJ - Warnings",
            anchor="cellranger-vdj-warnings",
            description="Warnings encountered during the analysis",
            plot=table.plot(
                warnings_by_sample,
                vdj_warnings_headers,
                {
                    "namespace": "VDJ",
                    "id": "cellranger-vdj-warnings-table",
                    "title": "Cellranger VDJ: Warnings",
                },
            ),
        )

    module.add_section(
        name="VDJ - Summary stats",
        anchor="cellranger-vdj-stats",
        description="Summary QC metrics from Cell Ranger count",
        plot=table.plot(
            mapping_by_sample,
            mapping_headers,
            {
                "namespace": "VDJ",
                "id": "cellranger-vdj-stats-table",
                "title": "Cellranger VDJ: Summary stats",
            },
        ),
    )

    module.add_section(
        name="VDJ - Annotations",
        anchor="cellranger-vdj-annot",
        description="V(D)J annotations from Cell Ranger VDJ analysis",
        plot=table.plot(
            annotations_by_sample,
            annotations_headers,
            {
                "namespace": "VDJ",
                "id": "cellranger-vdj-annot-table",
                "title": "Cellranger VDJ: Annotations",
            },
        ),
    )

    module.add_section(
        name="VDJ - BC rank plot",
        anchor="cellranger-vdj-bcrank-plot",
        description=plots_params_by_id["bc"]["description"],
        helptext=plots_params_by_id["bc"]["helptext"],
        plot=linegraph.plot(plots_data_by_sample_by_id["bc"], plots_params_by_id["bc"]["config"]),
    )

    return len(general_data_by_sample)
