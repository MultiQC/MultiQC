"""MultiQC module to parse output from Space Ranger count"""

import json
import logging
import os
from collections import defaultdict
from typing import Dict, Union

from multiqc import config, BaseMultiqcModule
from multiqc.plots import linegraph, table

from multiqc.modules.spaceranger.utils import set_hidden_cols, transform_data, populate_data_and_headers

log = logging.getLogger(__name__)


def parse_count_html(module: BaseMultiqcModule):
    """
    Space Ranger count report parser
    """

    general_stats_data: Dict[str, Dict[str, Union[str, float, int, None]]] = defaultdict()
    summary_data_by_sample: Dict[str, Dict[str, Union[str, float, int, None]]] = defaultdict()
    warnings_data_by_sample: Dict[str, Dict[str, Union[str, float, int, None]]] = defaultdict(lambda: defaultdict())

    plots_data: Dict[str, Dict] = {"saturation": defaultdict(), "genes": defaultdict(), "genomic_dna": defaultdict()}
    plots_conf: Dict[str, Dict] = {"saturation": defaultdict(), "genes": defaultdict(), "genomic_dna": defaultdict()}

    warnings_headers: Dict = dict()
    summary_headers: Dict[str, Dict[str, Union[str, float, int, None]]] = {
        "reads": {
            "rid": "count_data_reads",
            "title": "Reads",
            "description": "Number of reads",
            "shared_key": "read_count",
        }
    }
    general_stats_headers: Dict[str, Dict[str, Union[str, float, int, None]]] = {
        "reads": {
            "rid": "count_genstats_reads",
            "title": "Reads",
            "description": "Number of reads",
            "shared_key": "read_count",
            "namespace": "Space Ranger Count",
        }
    }

    for f in module.find_log_files("spaceranger/count_html", filehandles=True):
        # Go through the html report of space ranger and extract the data in a dicts
        summary = None
        for line in f["f"]:
            line = line.strip()
            if line.startswith("const data"):
                line = line.replace("const data = ", "")
                summary = json.loads(line)
                summary = summary["summary"]
                break

        if summary is None:
            logging.error(f"Couldn't find JSON summary data in HTML report, skipping: {f['fn']}")
            continue

        sample_name = module.clean_s_name(summary["sample"]["id"], f)
        if sample_name in general_stats_data:
            log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {sample_name}")

        software = next(
            iter(x[1] for x in summary["summary_tab"]["pipeline_info_table"]["rows"] if x[0] == "Pipeline Version")
        )
        software_name, software_version = software.split("-")
        module.add_software_version(version=software_version, sample=sample_name, software_name=software_name)

        # List of data collated from different tables in cellranger reports.
        # This is a list of Tuples (metric name, value)
        data_rows = (
            [["Number of Spots Under Tissue", summary["summary_tab"]["filtered_bcs_transcriptome_union"]["metric"]]]
            + summary["summary_tab"]["cells"]["table"]["rows"]
            + summary["summary_tab"]["sequencing"]["table"]["rows"]
            + summary["summary_tab"]["mapping"]["table"]["rows"]
        )
        # This is only contained in spaceranger reports for analyses with probe sets
        try:
            data_rows.extend(summary["analysis_tab"]["gdna"]["gems"]["table"]["rows"])
        except KeyError as e:
            fname = os.path.join(f["root"], f["fn"])
            log.debug(
                f"Could not parse the expected DNA table in the spaceranger report: {fname}: '{e}' field is missing"
            )

        general_stats_data[sample_name] = populate_data_and_headers(
            headers_to_update=general_stats_headers,
            new_data=data_rows,
            new_headers={
                "Number of Spots Under Tissue": "spots under tissue",
                "Mean Reads per Spot": "avg reads/spot",
                "Fraction Reads in Spots Under Tissue": "reads in spots",
                "Number of Reads": "reads",
                "Valid Barcodes": "valid bc",
            },
            colors={
                "spots under tissue": "RdPu",
                "avg reads/spot": "Blues",
                "reads in spots": "PiYG",
                "reads": "YlGn",
                "valid bc": "RdYlGn",
            },
            prefix="Count",
            int_cols=[
                "reads",
                "spots under tissue",
                "avg reads/spot",
            ],
        )
        if len(general_stats_data[sample_name]) == 0:
            continue

        summary_data_by_sample[sample_name] = populate_data_and_headers(
            headers_to_update=summary_headers,
            new_data=data_rows,
            new_headers={
                "Number of Reads": "reads",
                "Number of Spots Under Tissue": "spots under tissue",
                "Mean Reads per Spot": "avg reads/spot",
                "Fraction Reads in Spots Under Tissue": "reads in spots",
                "Genes Detected": "genes detected",
                "Median Genes per Spot": "median genes/spot",
                "Median UMI Counts per Spot": "median umi/spot",
                "Valid Barcodes": "valid bc",
                "Valid UMIs": "valid umi",
                "Sequencing Saturation": "saturation",
                "Q30 Bases in Barcode": "Q30 bc",
                "Q30 Bases in UMI": "Q30 UMI",
                "Q30 Bases in RNA Read": "Q30 read",
                "Reads Mapped to Probe Set": "reads mapped",
                "Reads Mapped Confidently to Probe Set": "confident reads",
                "Reads Mapped Confidently to Filtered Probe Set": "confident filtered reads",
                "Estimated UMIs from Genomic DNA": "genomic umis",
                "Estimated UMIs from Genomic DNA per Unspliced Probe": "genomic umis/unspliced probe",
            },
            colors={
                "reads": "YlGn",
                "spots under tissue": "RdPu",
                "avg reads/spot": "Blues",
                "genes detected": "Greens",
                "median genes/spot": "Purples",
                "reads in spots": "PuBuGn",
                "valid bc": "Spectral",
                "valid umi": "RdYlGn",
                "median umi/spot": "YlGn",
                "saturation": "PRGn",
                "genomic umis": "PuRd",
                "genomic umis/unspliced probe": "YlOrRd",
            },
            prefix="Count",
            int_cols=[
                "reads",
                "spots under tissue",
                "avg reads/spot",
                "median genes/spot",
                "genes detected",
                "genomic umis/unspliced probe",
            ],
        )
        if len(summary_data_by_sample[sample_name]) == 0:
            continue

        # Extract warnings if any
        alarms_list = summary["alarms"].get("alarms", [])
        for alarm in alarms_list:
            # "Intron mode used" alarm added in Space Ranger 7.0 lacks id
            if "id" not in alarm:
                continue
            warnings_data_by_sample[sample_name][alarm["id"]] = "FAIL"
            warnings_headers[alarm["id"]] = {
                "title": alarm["id"].replace("_", " ").title(),
                "description": alarm["title"],
                "bgcols": {"FAIL": "#f7dddc"},
            }

        # Extract data for plots
        # `analysis_tab` may not be present in the report if there are few reads
        try:
            if not plots_conf.get("saturation"):
                plots_conf["saturation"] = {
                    "config": {
                        "id": "mqc_spaceranger_count_saturation",
                        "title": f"Space Ranger: count: {summary['analysis_tab']['seq_saturation_plot']['help']['title']}",
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
            plots_data["saturation"][sample_name] = transform_data(
                summary["analysis_tab"]["seq_saturation_plot"]["plot"]["data"][0]
            )
        except KeyError as e:
            log.debug("No saturation plot found in the spaceranger report:", e)

        # `analysis_tab` may not be present in the report if there are few reads
        try:
            if not plots_conf.get("genes"):
                plots_conf["genes"] = {
                    "config": {
                        "id": "mqc_spaceranger_count_genesXspot",
                        "title": f"Space Ranger: count: {summary['analysis_tab']['median_gene_plot']['help']['title']}",
                        "xlab": summary["analysis_tab"]["median_gene_plot"]["plot"]["layout"]["xaxis"]["title"],
                        "ylab": summary["analysis_tab"]["median_gene_plot"]["plot"]["layout"]["yaxis"]["title"],
                        "ylog": False,
                        "xlog": False,
                    },
                    "description": "Median gene counts per spot",
                    "helptext": summary["analysis_tab"]["median_gene_plot"]["help"]["helpText"],
                }
            plots_data["genes"][sample_name] = transform_data(
                summary["analysis_tab"]["median_gene_plot"]["plot"]["data"][0]
            )
        except KeyError:
            log.debug("No median gene plot found in the spaceranger report")

        # The gDNA plots are only contained for spaceranger workflows with probesets
        try:
            if not plots_conf.get("genomic_dna"):
                plots_conf["genomic_dna"] = {
                    "config": {
                        "id": "mqc_spaceranger_count_genomic_dna",
                        "title": f"Space Ranger: count: {summary['analysis_tab']['gdna']['gems']['help']['title']}",
                        "xlab": summary["analysis_tab"]["gdna"]["plot"]["layout"]["xaxis"]["title"],
                        "ylab": summary["analysis_tab"]["gdna"]["plot"]["layout"]["yaxis"]["title"],
                        "ylog": False,
                        "xlog": False,
                    },
                    "description": "Estimated UMIs from Genomic DNA per Unspliced Probe",
                    "helptext": summary["analysis_tab"]["gdna"]["gems"]["help"]["data"][2][1][0]
                    + "\n\nThis summary graphic in the MultiQC report only shows the estimated mean baseline level of unspliced probe counts.",
                }
            plots_data["genomic_dna"][sample_name] = transform_data(summary["analysis_tab"]["gdna"]["plot"]["data"][2])
        except KeyError:
            log.debug("No genomic DNA plot found in the spaceranger report")

        module.add_data_source(f, sample_name, module="spaceranger", section="count")

    summary_data_by_sample = module.ignore_samples(summary_data_by_sample)
    warnings_data_by_sample = module.ignore_samples(warnings_data_by_sample)
    general_stats_data = module.ignore_samples(general_stats_data)
    for k in plots_data.keys():
        plots_data[k] = module.ignore_samples(plots_data[k])

    summary_headers = set_hidden_cols(
        summary_headers,
        [
            "Q30 bc",
            "Q30 UMI",
            "Q30 read",
            "reads in spots",
            "avg reads/spot",
            "confident reads",
            "confident filtered reads",
            "genomic unis/unspliced probe",
            "valid umi",
            "saturation",
        ],
    )

    module.general_stats_addcols(general_stats_data, general_stats_headers)

    # Write parsed report data to a file
    module.write_data_file(summary_data_by_sample, "multiqc_spaceranger_count")

    # Add sections to the report
    if len(warnings_data_by_sample) > 0:
        module.add_section(
            name="Count - Warnings",
            anchor="spaceranger-count-warnings-section",
            description="Warnings encountered during the analysis",
            plot=table.plot(
                warnings_data_by_sample,
                warnings_headers,
                {
                    "namespace": "Space Ranger Count",
                    "id": "spaceranger-count-warnings",
                    "title": "Space Ranger: Count: Warnings",
                },
            ),
        )

    module.add_section(
        name="Count - Summary stats",
        anchor="spaceranger-count-stats-section",
        description="Summary QC metrics from Space Ranger count",
        plot=table.plot(
            summary_data_by_sample,
            summary_headers,
            {
                "namespace": "Space Ranger Count",
                "id": "spaceranger-count-stats",
                "title": "Space Ranger: Count: Summary stats",
            },
        ),
    )

    # The gDNA plots are only contained for spaceranger workflows with probesets
    if plots_conf.get("genomic_dna"):
        module.add_section(
            name="Count - UMIs from Genomic DNA",
            anchor="spaceranger-count-bcrank-plot-section",
            description=plots_conf["genomic_dna"]["description"],
            helptext=plots_conf["genomic_dna"]["helptext"],
            plot=linegraph.plot(
                plots_data["genomic_dna"],
                plots_conf["genomic_dna"]["config"],
            ),
        )

    if plots_conf.get("genes"):
        module.add_section(
            name="Count - Median genes",
            anchor="spaceranger-count-genes-plot-section",
            description=plots_conf["genes"]["description"],
            helptext=plots_conf["genes"]["helptext"],
            plot=linegraph.plot(plots_data["genes"], plots_conf["genes"]["config"]),
        )

    if plots_conf.get("saturation"):
        module.add_section(
            name="Count - Saturation plot",
            anchor="spaceranger-count-saturation-plot-section",
            description=plots_conf["saturation"]["description"],
            helptext=plots_conf["saturation"]["helptext"],
            plot=linegraph.plot(
                plots_data["saturation"],
                plots_conf["saturation"]["config"],
            ),
        )

    return len(general_stats_data)
