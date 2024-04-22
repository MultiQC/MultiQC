""" MultiQC module to parse output from Space Ranger count """

import json
import logging
import os

from multiqc import config
from multiqc.plots import linegraph, table

from ._utils import set_hidden_cols, transform_data, update_data_and_headers

# Initialise the logger
log = logging.getLogger(__name__)


class SpaceRangerCountMixin:
    """Space Ranger count report parser"""

    def parse_count_html(self):
        self.generalstats_data = dict()
        self.summary_data = dict()
        self.warnings_data = dict()
        self.generalstats_headers = dict()
        self.summary_headers = dict()
        self.warnings_headers = dict()
        self.plots_conf = {"saturation": dict(), "genes": dict(), "genomic_dna": dict()}
        self.plots_data = {"saturation": dict(), "genes": dict(), "genomic_dna": dict()}

        for f in self.find_log_files("spaceranger/count_html", filehandles=True):
            self.parse_count_report(f)

        self.summary_data = self.ignore_samples(self.summary_data)
        self.generalstats_data = self.ignore_samples(self.generalstats_data)
        self.warnings_data = self.ignore_samples(self.warnings_data)
        for k in self.plots_data.keys():
            self.plots_data[k] = self.ignore_samples(self.plots_data[k])

        self.generalstats_headers["reads"] = {
            "rid": "count_genstats_reads",
            "title": "Reads",
            "description": f"Number of reads ({config.read_count_desc})",
            "shared_key": "read_count",
            "namespace": "Space Ranger Count",
        }

        self.summary_headers["reads"] = {
            "rid": "count_data_reads",
            "title": "Reads",
            "description": f"Number of reads ({config.read_count_desc})",
            "shared_key": "read_count",
        }
        self.summary_headers = set_hidden_cols(
            self.summary_headers,
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

        if len(self.generalstats_data) == 0:
            return 0

        self.general_stats_addcols(self.generalstats_data, self.generalstats_headers)

        # Write parsed report data to a file
        self.write_data_file(self.summary_data, "multiqc_spaceranger_count")

        # Add sections to the report
        if len(self.warnings_data) > 0:
            self.add_section(
                name="Count - Warnings",
                anchor="spaceranger-count-warnings-section",
                description="Warnings encountered during the analysis",
                plot=table.plot(
                    self.warnings_data,
                    self.warnings_headers,
                    {
                        "namespace": "Space Ranger Count",
                        "id": "spaceranger-count-warnings",
                        "title": "Space Ranger Count: Warnings",
                    },
                ),
            )

        self.add_section(
            name="Count - Summary stats",
            anchor="spaceranger-count-stats-section",
            description="Summary QC metrics from Space Ranger count",
            plot=table.plot(
                self.summary_data,
                self.summary_headers,
                {
                    "namespace": "Space Ranger Count",
                    "id": "spaceranger-count-stats",
                    "title": "Space Ranger Count: Summary stats",
                },
            ),
        )

        # The gDNA plots are only contained for spaceranger workflows with probesets
        if "genomic_dna" in self.plots_conf:
            self.add_section(
                name="Count - UMIs from Genomic DNA",
                anchor="spaceranger-count-bcrank-plot-section",
                description=self.plots_conf["genomic_dna"]["description"],
                helptext=self.plots_conf["genomic_dna"]["helptext"],
                plot=linegraph.plot(
                    self.plots_data["genomic_dna"],
                    self.plots_conf["genomic_dna"]["config"],
                ),
            )

        if "genes" in self.plots_conf:
            self.add_section(
                name="Count - Median genes",
                anchor="spaceranger-count-genes-plot-section",
                description=self.plots_conf["genes"]["description"],
                helptext=self.plots_conf["genes"]["helptext"],
                plot=linegraph.plot(self.plots_data["genes"], self.plots_conf["genes"]["config"]),
            )

        if "saturation" in self.plots_conf:
            self.add_section(
                name="Count - Saturation plot",
                anchor="spaceranger-count-saturation-plot-section",
                description=self.plots_conf["saturation"]["description"],
                helptext=self.plots_conf["saturation"]["helptext"],
                plot=linegraph.plot(
                    self.plots_data["saturation"],
                    self.plots_conf["saturation"]["config"],
                ),
            )

        return len(self.generalstats_data)

    def parse_count_report(self, f):
        """Go through the html report of space ranger and extract the data in a dicts"""

        summary = None
        for line in f["f"]:
            line = line.strip()
            if line.startswith("const data"):
                line = line.replace("const data = ", "")
                summary = json.loads(line)
                summary = summary["summary"]
                break

        assert summary is not None, "Couldn't find JSON summary data in HTML report."

        sample_name = self.clean_s_name(summary["sample"]["id"], f)
        software = next(
            iter(x[1] for x in summary["summary_tab"]["pipeline_info_table"]["rows"] if x[0] == "Pipeline Version")
        )
        software_name, software_version = software.split("-")
        self.add_software_version(version=software_version, sample=sample_name, software_name=software_name)

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

        # Store general stats
        col_dict = {
            "Number of Spots Under Tissue": "spots under tissue",
            "Mean Reads per Spot": "avg reads/spot",
            "Fraction Reads in Spots Under Tissue": "reads in spots",
            "Number of Reads": "reads",
            "Valid Barcodes": "valid bc",
        }
        colors = {
            "spots under tissue": "RdPu",
            "avg reads/spot": "Blues",
            "reads in spots": "PiYG",
            "reads": "YlGn",
            "valid bc": "RdYlGn",
        }
        int_cols = [
            "reads",
            "spots under tissue",
            "avg reads/spot",
        ]
        data_general_stats = {}
        update_data_and_headers(
            data=data_general_stats,
            headers=self.generalstats_headers,
            new_data=data_rows,
            new_headers=col_dict,
            colors=colors,
            prefix="Count",
            int_cols=int_cols,
        )

        # Store full data from space ranger count report
        col_dict = {
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
        }
        colors = {
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
        }
        int_cols = [
            "reads",
            "spots under tissue",
            "avg reads/spot",
            "median genes/spot",
            "genes detected",
            "genomic umis/unspliced probe",
        ]
        data = {}
        update_data_and_headers(
            data=data,
            headers=self.summary_headers,
            new_data=data_rows,
            new_headers=col_dict,
            colors=colors,
            prefix="Count",
            int_cols=int_cols,
        )

        # Extract warnings if any
        warnings = {}
        alarms_list = summary["alarms"].get("alarms", [])
        for alarm in alarms_list:
            # "Intron mode used" alarm added in Space Ranger 7.0 lacks id
            if "id" not in alarm:
                continue
            warnings[alarm["id"]] = "FAIL"
            self.warnings_headers[alarm["id"]] = {
                "title": alarm["id"].replace("_", " ").title(),
                "description": alarm["title"],
                "bgcols": {"FAIL": "#f7dddc"},
            }

        # Extract data for plots
        plots_conf = {}
        plots_data = {}
        # `analysis_tab` may not be present in the report if there are few reads
        try:
            plots_conf["saturation"] = {
                "config": {
                    "id": "mqc_spaceranger_count_saturation",
                    "title": f"Space Ranger count: {summary['analysis_tab']['seq_saturation_plot']['help']['title']}",
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
            plots_data["saturation"] = {
                sample_name: transform_data(summary["analysis_tab"]["seq_saturation_plot"]["plot"]["data"][0])
            }
        except KeyError as e:
            log.debug("No saturation plot found in the spaceranger report:", e)
            if "saturation" in plots_conf:
                del plots_conf["saturation"]
            if "saturation" in plots_data:
                del plots_data["saturation"]

        # `analysis_tab` may not be present in the report if there are few reads
        try:
            plots_conf["genes"] = {
                "config": {
                    "id": "mqc_spaceranger_count_genesXspot",
                    "title": f"Space Ranger count: {summary['analysis_tab']['median_gene_plot']['help']['title']}",
                    "xlab": summary["analysis_tab"]["median_gene_plot"]["plot"]["layout"]["xaxis"]["title"],
                    "ylab": summary["analysis_tab"]["median_gene_plot"]["plot"]["layout"]["yaxis"]["title"],
                    "yLog": False,
                    "xLog": False,
                },
                "description": "Median gene counts per spot",
                "helptext": summary["analysis_tab"]["median_gene_plot"]["help"]["helpText"],
            }
            plots_data["genes"] = {
                sample_name: transform_data(summary["analysis_tab"]["median_gene_plot"]["plot"]["data"][0])
            }
        except KeyError:
            log.debug("No median gene plot found in the spaceranger report")
            if "genes" in plots_conf:
                del plots_conf["genes"]
            if "genes" in plots_data:
                del plots_data["genes"]

        # The gDNA plots are only contained for spaceranger workflows with probesets
        try:
            plots_conf["genomic_dna"] = {
                "config": {
                    "id": "mqc_spaceranger_count_genomic_dna",
                    "title": f"Space Ranger count: {summary['analysis_tab']['gdna']['gems']['help']['title']}",
                    "xlab": summary["analysis_tab"]["gdna"]["plot"]["layout"]["xaxis"]["title"],
                    "ylab": summary["analysis_tab"]["gdna"]["plot"]["layout"]["yaxis"]["title"],
                    "yLog": False,
                    "xLog": False,
                },
                "description": "Estimated UMIs from Genomic DNA per Unspliced Probe",
                "helptext": summary["analysis_tab"]["gdna"]["gems"]["help"]["data"][2][1][0]
                + "\n\nThis summary graphic in the MultiQC report only shows the estimated mean baseline level of unspliced probe counts.",
            }
            plots_data["genomic_dna"] = {
                sample_name: transform_data(summary["analysis_tab"]["gdna"]["plot"]["data"][2])
            }
        except KeyError:
            log.debug("No genomic DNA plot found in the spaceranger report")
            if "genomic_dna" in plots_conf:
                del plots_conf["genomic_dna"]
            if "genomic_dna" in plots_data:
                del plots_data["genomic_dna"]

        if len(data) > 0:
            if sample_name in self.generalstats_data:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {sample_name}")
            self.add_data_source(f, sample_name, module="spaceranger", section="count")
            self.summary_data[sample_name] = data
            self.generalstats_data[sample_name] = data_general_stats
            if len(warnings) > 0:
                self.warnings_data[sample_name] = warnings
            self.plots_conf = plots_conf
            for k in plots_data.keys():
                if k not in self.plots_data.keys():
                    self.plots_data[k] = dict()
                self.plots_data[k].update(plots_data[k])
