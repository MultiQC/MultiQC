""" MultiQC module to parse output from Space Ranger count """

import json
import logging
from collections import OrderedDict

from multiqc import config
from multiqc.plots import linegraph, table

from ._utils import set_hidden_cols, transform_data, update_dict

# Initialise the logger
log = logging.getLogger(__name__)


class SpaceRangerCountMixin:
    """Space Ranger count report parser"""

    def parse_count_html(self):
        self.spacerangercount_data = dict()
        self.spacerangercount_general_data = dict()
        self.spacerangercount_warnings = dict()
        self.spacerangercount_plots_conf = {"saturation": dict(), "genes": dict(), "genomic_dna": dict()}
        self.spacerangercount_plots_data = {"saturation": dict(), "genes": dict(), "genomic_dna": dict()}
        self.count_general_data_headers = dict()
        self.count_data_headers = OrderedDict()
        self.count_warnings_headers = OrderedDict()

        for f in self.find_log_files("spaceranger/count_html", filehandles=True):
            self.parse_count_report(f)

        self.spacerangercount_data = self.ignore_samples(self.spacerangercount_data)
        self.spacerangercount_general_data = self.ignore_samples(self.spacerangercount_general_data)
        self.spacerangercount_warnings = self.ignore_samples(self.spacerangercount_warnings)
        for k in self.spacerangercount_plots_data.keys():
            self.spacerangercount_plots_data[k] = self.ignore_samples(self.spacerangercount_plots_data[k])

        self.count_general_data_headers["reads"] = {
            "rid": "count_genstats_reads",
            "title": "{} Reads".format(config.read_count_prefix),
            "description": "Number of reads ({})".format(config.read_count_desc),
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "namespace": "Space Ranger Count",
        }

        self.count_data_headers["reads"] = {
            "rid": "count_data_reads",
            "title": "{} Reads".format(config.read_count_prefix),
            "description": "Number of reads ({})".format(config.read_count_desc),
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.count_data_headers = set_hidden_cols(
            self.count_data_headers,
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

        if len(self.spacerangercount_general_data) == 0:
            return 0

        else:
            self.general_stats_addcols(self.spacerangercount_general_data, self.count_general_data_headers)

            # Write parsed report data to a file
            self.write_data_file(self.spacerangercount_data, "multiqc_spaceranger_count")

            # Add sections to the report
            if len(self.spacerangercount_warnings) > 0:
                self.add_section(
                    name="Count - Warnings",
                    anchor="spaceranger-count-warnings",
                    description="Warnings encountered during the analysis",
                    plot=table.plot(
                        self.spacerangercount_warnings,
                        self.count_warnings_headers,
                        {"namespace": "Space Ranger Count", "id": "spaceranger-count-warnings"},
                    ),
                )

            self.add_section(
                name="Count - Summary stats",
                anchor="spaceranger-count-stats",
                description="Summary QC metrics from Space Ranger count",
                plot=table.plot(
                    self.spacerangercount_data,
                    self.count_data_headers,
                    {"namespace": "Space Ranger Count", "id": "spaceranger-count-stats"},
                ),
            )

            # The gDNA plots are only contained for spaceranger workflows with probesets
            try:
                self.add_section(
                    name="Count - UMIs from Genomic DNA",
                    anchor="spaceranger-count-bcrank-plot",
                    description=self.spacerangercount_plots_conf["genomic_dna"]["description"],
                    helptext=self.spacerangercount_plots_conf["genomic_dna"]["helptext"],
                    plot=linegraph.plot(
                        self.spacerangercount_plots_data["genomic_dna"],
                        self.spacerangercount_plots_conf["genomic_dna"]["config"],
                    ),
                )
            except KeyError:
                pass

            self.add_section(
                name="Count - Median genes",
                anchor="spaceranger-count-genes-plot",
                description=self.spacerangercount_plots_conf["genes"]["description"],
                helptext=self.spacerangercount_plots_conf["genes"]["helptext"],
                plot=linegraph.plot(
                    self.spacerangercount_plots_data["genes"], self.spacerangercount_plots_conf["genes"]["config"]
                ),
            )

            self.add_section(
                name="Count - Saturation plot",
                anchor="spaceranger-count-saturation-plot",
                description=self.spacerangercount_plots_conf["saturation"]["description"],
                helptext=self.spacerangercount_plots_conf["saturation"]["helptext"],
                plot=linegraph.plot(
                    self.spacerangercount_plots_data["saturation"],
                    self.spacerangercount_plots_conf["saturation"]["config"],
                ),
            )

            return len(self.spacerangercount_general_data)

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
        # This is only contained in spaceranger reports for anlayses with probesets
        try:
            data_rows.extend(summary["analysis_tab"]["gdna"]["gems"]["table"]["rows"])
        except KeyError:
            pass

        # Store general stats
        col_dict = {
            "Number of Spots Under Tissue": "spots under tissue",
            "Mean Reads per Spot": "avg reads/spot",
            "Fraction Reads in Spots Under Tissue": "reads in spots",
            "Number of Reads": "reads",
            "Valid Barcodes": "valid bc",
        }
        colours = {
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
        update_dict(
            data_general_stats, self.count_general_data_headers, data_rows, col_dict, colours, "Count", int_cols
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
        colours = {
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
        update_dict(data, self.count_data_headers, data_rows, col_dict, colours, "Count", int_cols)

        # Extract warnings if any
        warnings = {}
        alarms_list = summary["alarms"].get("alarms", [])
        for alarm in alarms_list:
            # "Intron mode used" alarm added in Space Ranger 7.0 lacks id
            if "id" not in alarm:
                continue
            warnings[alarm["id"]] = "FAIL"
            self.count_warnings_headers[alarm["id"]] = {
                "title": alarm["id"],
                "description": alarm["title"],
                "bgcols": {"FAIL": "#f06807"},
            }

        # Extract data for plots
        plots = {}
        plots_data = {}
        # `analysis_tab` may not be present in the report if there are few reads
        try:
            plots["saturation"] = {
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
        except KeyError:
            pass

        # `analysis_tab` may not be present in the report if there are few reads
        try:
            plots["genes"] = {
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
            pass

        # The gDNA plots are only contained for spaceranger workflows with probesets
        try:
            plots["genomic_dna"] = {
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
            pass

        if len(data) > 0:
            if sample_name in self.spacerangercount_general_data:
                log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], sample_name))
            self.add_data_source(f, sample_name, module="spaceranger", section="count")
            self.spacerangercount_data[sample_name] = data
            self.spacerangercount_general_data[sample_name] = data_general_stats
            if len(warnings) > 0:
                self.spacerangercount_warnings[sample_name] = warnings
            self.spacerangercount_plots_conf = plots
            for k in plots_data.keys():
                if k not in self.spacerangercount_plots_data.keys():
                    self.spacerangercount_plots_data[k] = dict()
                self.spacerangercount_plots_data[k].update(plots_data[k])
