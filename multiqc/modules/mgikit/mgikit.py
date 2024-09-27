#!/usr/bin/env python

"""MultiQC example plugin module"""

from __future__ import print_function
from collections import OrderedDict
from multiqc import config
from multiqc.plots import table, bargraph
from multiqc.base_module import BaseMultiqcModule
#from multiqc.base_module import ModuleNoSamplesFound
import logging

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name="MGIKIT",
            target="mgikit",
            anchor="mgikit",
            href="https://github.com/sagc-bioinformatics/mgikit",
            info="can be used to demultiplex demultiplex FASTQ files from an MGI sequencing instrument for downstream analysis.",
            doi="https://doi.org/10.1101/2024.01.09.574938",
        )

        # Halt execution if we've disabled the plugin
        if config.kwargs.get("mgikit_disable_plugin", True):
            return None

        for f in self.find_log_files("mgikit"):
            self.add_data_source(f)

        self.add_software_version(None)
        module_name = "MGIKIT"
        # Colours from Set3 of ColorBrewer R package, with a few additions for perfect matches (#1f78b4),
        # undetermined reads (#000000) and ambiguous reads (#a6761d)
        mgikit_colors = [
            "#1f78b4",
            "#8dd3c7",
            "#ffffb3",
            "#bebada",
            "#fb8072",
            "#80b1d3",
            "#fdb462",
            "#b3de69",
            "#fccde5",
            "#d9d9d9",
            "#bc80bd",
            "#ccebc5",
            "#ffed6f",
            "#000000",
            "#a6761d",
        ]
        undetermined_label = "Undetermined"
        ambiguous_label = "Ambiguous"

        visualisation_threshold = config.kwargs.get("undetermined_barcode_threshold", 25)
        brief_report = config.kwargs.get("brief_report", False)
        if brief_report:
            log.warning("The brief version version is activated!")

        show_all_samples = not config.kwargs.get("keep_core_samples", False)
        if not show_all_samples:
            log.warning("Undetermined and Ambiguous cases will be eliminated in this report!")
        # Find and load any input files for this module

        decimal_positions = config.kwargs.get("decimal_positions", 2)
        decimal_positions = str(decimal_positions)
        ################
        # General statistics
        general_info_logs = self.find_log_files("mgikit/mgi_general_info")
        mgi_general_statistics = {}
        columns_headers = OrderedDict()
        columns_headers["M Clusters"] = {
            "namespace": "mgikit",
            "title": "M Total Clusters",  # Short title, table column title
            "description": "Total number of clusters for this sample (millions)",
            "format": "{:,." + decimal_positions + "f}",
        }
        columns_headers["Mb Yield"] = {
            "namespace": "mgikit",
            "title": "Mb Yield",  # Short title, table column title
            "description": "Number of bases (millions)",
            "format": "{:,." + decimal_positions + "f}",
        }
        columns_headers["Mb Yield ≥ Q30"] = {
            "namespace": "mgikit",
            "title": "Mb Yield ≥ Q30",  # Short title, table column title
            "description": "Number of bases (millions)",
            "format": "{:,." + decimal_positions + "f}",
        }
        columns_headers["% R1 Yield ≥ Q30"] = {
            "namespace": "mgikit",
            "title": "% R1 Yield ≥ Q30",  # Short title, table column title
            "description": "Percent of bases in R1 with phred qualty score  ≥ Q30",
            "max": 100,  # Minimum value in range, for bar / colour coding
            "min": 0,  # Maximum value in range, for bar / colour coding
            "suffix": "%",  # Suffix for value (eg. '%')
            "format": "{:,." + decimal_positions + "f}",
        }
        columns_headers["% R2 Yield ≥ Q30"] = {
            "namespace": "mgikit",
            "title": "% R2 Yield ≥ Q30",  # Short title, table column title
            "description": "Percent of bases in R2 with phred qualty score  ≥ Q30",
            "max": 100,  # Minimum value in range, for bar / colour coding
            "min": 0,  # Maximum value in range, for bar / colour coding
            "suffix": "%",  # Suffix for value (eg. '%')
            "format": "{:,." + decimal_positions + "f}",
        }
        columns_headers["% R3 Yield ≥ Q30"] = {
            "namespace": "mgikit",
            "title": "% R3 Yield ≥ Q30",  # Short title, table column title
            "description": "Percent of bases in R3 with phred qualty score  ≥ Q30",
            "max": 100,  # Minimum value in range, for bar / colour coding
            "min": 0,  # Maximum value in range, for bar / colour coding
            "suffix": "%",  # Suffix for value (eg. '%')
            "hidden": True,
            "format": "{:,." + decimal_positions + "f}",
        }
        columns_headers["% Perfect Index"] = {
            "namespace": "mgikit",
            "title": "% Perfect Index",  # Short title, table column title
            "description": "Percent of reads with perfect index (0 mismatches)",
            # Longer description, goes in mouse hover text
            "max": 100,  # Minimum value in range, for bar / colour coding
            "min": 0,  # Maximum value in range, for bar / colour coding
            "suffix": "%",  # Suffix for value (eg. '%')
            "format": "{:,." + decimal_positions + "f}",
        }
        file_cnt = 0
        for f in general_info_logs:
            file_cnt += 1
            collect_data = False
            lines = f["f"].splitlines()
            line_itr = 0
            header = None
            while line_itr < len(lines):
                line = lines[line_itr]
                if collect_data:
                    if len(line) == 0:
                        break
                    else:
                        vals = line.split("\t")
                        mgi_general_statistics[vals[0]] = {}
                        for i in range(1, len(header)):
                            mgi_general_statistics[vals[0]][header[i]] = vals[i]

                if line.startswith("#sample general info"):
                    collect_data = True
                    line_itr += 1
                    header = lines[line_itr].split("\t")
                line_itr += 1

        # print(mgi_general_statistics)
        if len(mgi_general_statistics.keys()) > 0:
            self.general_stats_addcols(mgi_general_statistics, columns_headers)
            self.write_data_file(mgi_general_statistics, "multiqc_mgikit_general")

        if file_cnt > 0:
            log.info("{} general information log files (*.mgikit.general) were loaded!".format(file_cnt))

        # lane statistics
        general_info_logs = self.find_log_files("mgikit/mgi_general_info")
        mgi_lane_statistics = {}
        columns_headers = OrderedDict()

        columns_headers["Mb Total Yield"] = {
            "namespace": "mgikit",
            "title": "Mb Total Yield",  # Short title, table column title
            "description": "Number of bases (millions)",
            # Longer description, goes in mouse hover text
            "min": 0,  # Maximum value in range, for bar / colour coding
            "format": "{:,." + decimal_positions + "f}",
        }

        columns_headers["M Total Clusters"] = {
            "namespace": "mgikit",
            "title": "M Total Clusters",  # Short title, table column title
            "description": "Total number of clusters for this lane (millions)",
            # Longer description, goes in mouse hover text
            "min": 0,  # Maximum value in range, for bar / colour coding
            "format": "{:,." + decimal_positions + "f}",
        }

        columns_headers["% bases ≥ Q30"] = {
            "namespace": "mgikit",
            "title": "% bases ≥ Q30",  # Short title, table column title
            "description": "Number of bases (millions)",
            # Longer description, goes in mouse hover text
            "max": 100,  # Minimum value in range, for bar / colour coding
            "min": 0,  # Maximum value in range, for bar / colour coding
            "suffix": "%",  # Suffix for value (eg. '%')
            "format": "{:,." + decimal_positions + "f}",
        }

        columns_headers["Mean Quality"] = {
            "namespace": "mgikit",
            "title": "Mean Quality",  # Short title, table column title
            "description": "Average phred qualty score",
            # Longer description, goes in mouse hover text
            "max": 40,  # Minimum value in range, for bar / colour coding
            "min": 0,  # Maximum value in range, for bar / colour coding
            "format": "{:,." + decimal_positions + "f}",
        }

        columns_headers["% Perfect Index"] = {
            "namespace": "mgikit",
            "title": "% Perfect Index",  # Short title, table column title
            "description": "Percent of reads with perfect index (0 mismatches)",
            # Longer description, goes in mouse hover text
            "max": 100,  # Minimum value in range, for bar / colour coding
            "min": 0,  # Maximum value in range, for bar / colour coding
            "suffix": "%",  # Suffix for value (eg. '%')
            "format": "{:,." + decimal_positions + "f}",
        }

        for f in general_info_logs:
            lines = f["f"].splitlines()
            for line_itr in range(len(lines)):
                line = lines[line_itr]
                if line.startswith("#Lane statistics"):
                    # found it. Information is in the next 2 lines
                    header = lines[line_itr + 1].split("\t")
                    vals = lines[line_itr + 2].split("\t")
                    mgi_lane_statistics[vals[0]] = {}
                    for i in range(1, len(header)):
                        mgi_lane_statistics[vals[0]][header[i]] = vals[i]
                    break
        if len(mgi_lane_statistics.keys()) > 0:
            self.add_section(
                name="Lane Statistics",
                anchor="mgikit-lane-statistics",
                description="Statistics about each lane for each flowcell.",
                helptext="This longer string (can be **markdown**) helps explain how to interpret the plot",
                plot=table.plot(
                    mgi_lane_statistics,
                    columns_headers,
                    {
                        "id": "mgi_lane_stats_plot",
                        "title": "Lane Statistics",
                        "sort_rows": False,
                        "col1_header": "Run ID - Lane",
                    },
                ),
            )

        # statistics about matching reads
        mgi_sample_reads_data = dict()
        mgi_lane_read_data = {}
        mgi_lane_sample_read_data = {}
        sample_read_logs = self.find_log_files("mgikit/mgi_sample_reads")
        file_cnt = 0
        for f in sample_read_logs:
            file_cnt += 1
            log_details = f["fn"].split(".")
            if len(log_details[-4]) == 0:
                curr_label = log_details[-3]
            else:
                curr_label = log_details[-4] + "-" + log_details[-3]

            # print(log_details)
            mgi_lane_read_data[curr_label] = {}
            mgi_lane_sample_read_data[curr_label] = {}
            lines = f["f"].splitlines()
            header = [x.strip() for x in lines[0].split("\t")]
            for i in range(1, len(header)):
                mgi_lane_read_data[curr_label][header[i]] = 0

            for line in lines[1:]:
                vals = line.split("\t")
                if len(vals) < 2 or vals[0] == "Total":
                    continue

                if vals[0] in [undetermined_label, ambiguous_label]:
                    mgi_lane_read_data[curr_label][vals[0]] = int(vals[1])

                if not show_all_samples and vals[0] in [undetermined_label, ambiguous_label]:
                    continue
                # print(vals[0], mgi_lane_read_data)
                if vals[0] in mgi_lane_sample_read_data[curr_label].keys():
                    log.warning("Same sample ({}) appears twice in the same log file!".format(vals[0]))

                mgi_lane_sample_read_data[curr_label][vals[0]] = {}

                if vals[0] not in mgi_sample_reads_data.keys():
                    mgi_sample_reads_data[vals[0]] = {}

                for i in range(1, len(header)):
                    if header[i] not in mgi_sample_reads_data[vals[0]].keys():
                        mgi_sample_reads_data[vals[0]][header[i]] = 0

                    mgi_lane_sample_read_data[curr_label][vals[0]][header[i]] = 0
                    mgi_sample_reads_data[vals[0]][header[i]] += int(vals[i])
                    if vals[0] not in [undetermined_label, ambiguous_label]:
                        mgi_lane_read_data[curr_label][header[i]] += int(vals[i])
                    mgi_lane_sample_read_data[curr_label][vals[0]][header[i]] += int(vals[i])
                # print(mgi_sample_reads_data)

        mgi_sample_reads_data_filtered = None
        if file_cnt > 0:
            log.info("{} Read information log files (*.mgikit.info) were loaded!".format(file_cnt))

            # Filter out samples matching ignored sample names
            mgi_sample_reads_data_filtered = self.ignore_samples(mgi_sample_reads_data)
            cats = OrderedDict()
            for i in range(1, len(header)):
                cats[header[i]] = {"name": header[i], "color": mgikit_colors[i - 1]}

            cats[undetermined_label] = {"name": undetermined_label, "color": mgikit_colors[len(header) - 1]}
            cats[ambiguous_label] = {"name": ambiguous_label, "color": mgikit_colors[len(header)]}
            # print(cats)
            # print(mgi_lane_read_data)
            pconfig = {
                "id": "mgi_lane_read_plot",
                "title": module_name + ": Clusters by lane",
                "xlab": "Lane",
                "ylab": "# Reads",
                "hide_zero_cats": True,
            }

            # print(mgi_lane_read_data)
            # print(cats)
            lane_read_plt = bargraph.plot(mgi_lane_read_data, cats, pconfig=pconfig)

            # Add a report section with the line plot
            self.add_section(
                name="Clusters by lane",
                anchor="mgikit-first",
                description="Number of reads per lane.",
                helptext="The number of reads with different amounts of mismatches in their indexes, up to the maximum allowed during demultiplexing.",
                plot=lane_read_plt,
            )

            # Nothing found - raise a UserWarning to tell MultiQC
            if not mgi_sample_reads_data_filtered or len(mgi_sample_reads_data_filtered.keys()) == 0:
                log.info("No sample left after filteration!")
                raise UserWarning
            else:
                pconfig = {
                    "id": "mgi_sample_read_plot",
                    "title": module_name + ": Clusters by sample",
                    "xlab": "Sample",
                    "ylab": "Matching reads",
                    "hide_zero_cats": True,
                }

                sample_read_plt = bargraph.plot(mgi_sample_reads_data_filtered, cats, pconfig=pconfig)

                # Add a report section with the line plot
                self.add_section(
                    name="Clusters by sample",
                    anchor="mgikit-second",
                    description="Number of reads per sample.",
                    helptext="Perfect index reads are those that do not have a single mismatch. All samples are aggregated across lanes combinned. Undetermined reads are treated as a separate sample.",
                    plot=sample_read_plt,
                )

                if not brief_report:
                    cntr = 0
                    for lane in sorted(mgi_lane_read_data.keys()):
                        cntr += 1
                        pconfig = {
                            "id": "mgi_sample_read_plot" + str(cntr),
                            "title": module_name + ": Clusters by sample for Lane (" + lane + ")",
                            "xlab": "Sample",
                            "ylab": "Matching reads",
                            "hide_zero_cats": True,
                        }

                        sample_read_plt = bargraph.plot(mgi_lane_sample_read_data[lane], cats, pconfig=pconfig)

                        # Add a report section with the line plot
                        self.add_section(
                            name="Clusters by sample for Lane: " + lane,
                            anchor="mgikit-lane-sample_reads-" + str(cntr),
                            description="This plot shows the number of reads per sample for Lane: " + lane,
                            helptext="This longer string (can be **markdown**) helps explain how to interpret the plot",
                            plot=sample_read_plt,
                        )

        # Check undetermined barcodes
        undetermined_barcode_logs = self.find_log_files("mgikit/mgi_undetermined_barcode")
        mgi_undetermined_barcode_data = {}
        mgi_undetermined_barcode_cnt = {}

        file_cnt = 0
        cat_lane = {}
        for f in undetermined_barcode_logs:
            # print(f['fn'])
            file_cnt += 1
            log_details = f["fn"].split(".")
            if len(log_details[-4]) == 0:
                curr_label = log_details[-3]
            else:
                curr_label = log_details[-4] + "-" + log_details[-3]

            lines = f["f"].splitlines()
            cat_lane[curr_label] = {
                "name": curr_label,
                "color": mgikit_colors[file_cnt - 1],
            }
            for line in lines:
                vals = [x.strip() for x in line.split()]
                if len(vals) < 2:
                    continue
                if vals[0] not in mgi_undetermined_barcode_data.keys():
                    mgi_undetermined_barcode_data[vals[0]] = {}
                mgi_undetermined_barcode_data[vals[0]][curr_label] = int(vals[1])

                if vals[0] not in mgi_undetermined_barcode_cnt.keys():
                    mgi_undetermined_barcode_cnt[vals[0]] = int(vals[1])
                else:
                    mgi_undetermined_barcode_cnt[vals[0]] += int(vals[1])

        if file_cnt > 0:
            log.info("{} Undetermined barcode log files (*.mgikit.undetermined_barcode) were loaded!".format(file_cnt))

            pconfig = {
                "id": "mgi_lane_undetermined_plot",
                "title": module_name + ": Undetermined barcodes by lane",
                "xlab": "Barcode",
                "ylab": "# Reads",
                "hide_zero_cats": True,
            }

            sorted_dic = OrderedDict(
                [
                    (k[0], mgi_undetermined_barcode_data[k[0]])
                    for k in sorted(
                        mgi_undetermined_barcode_cnt.items(),
                        key=lambda item: -int(item[1]),
                    )[:visualisation_threshold]
                ]
            )

            lane_barcode_plt = bargraph.plot(sorted_dic, cat_lane, pconfig=pconfig)

            # Add a report section with the line plot
            self.add_section(
                name="Undetermined barcodes by lane",
                anchor="mgikit-undetermined-barcode",
                description="Count of the top twenty five most abundant undetermined barcodes by lanes.",
                plot=lane_barcode_plt,
            )
        # else:
        #    log.warning("No log files (*.mgikit.undetermined_barcode) was found!")

        # Check ambiguous barcodes
        ambiguous_barcode_logs = self.find_log_files("mgikit/mgi_ambiguous_barcode")
        mgi_ambiguous_barcode_data = {}
        file_cnt = 0
        cat_lane = {}
        for f in ambiguous_barcode_logs:
            # print(f['fn'])
            file_cnt += 1
            log_details = f["fn"].split(".")
            if len(log_details[-4]) == 0:
                curr_label = log_details[-3]
            else:
                curr_label = log_details[-4] + "-" + log_details[-3]
            lines = f["f"].splitlines()
            for line in lines:
                vals = [x.strip() for x in line.split()]
                if len(vals) < 2:
                    continue
                mgi_ambiguous_barcode_data[curr_label + ":" + vals[0]] = {
                    "Frequency": int(vals[1]),
                    # "Samples": vals[2],
                }

        if file_cnt > 0:
            log.info("{} Ambiguous barcode log files (*.mgikit.ambiguous_barcode) were loaded!".format(file_cnt))
            sorted_dic = OrderedDict(
                [
                    (k, v)
                    for k, v in sorted(
                        mgi_ambiguous_barcode_data.items(),
                        key=lambda item: -int(item[1]["Frequency"]),
                    )
                ][:visualisation_threshold]
            )

            # Add a report section with the line plot
            columns_headers = {}
            columns_headers["Frequency"] = {
                "namespace": "mgikit",
                "title": "Frequency",  # Short title, table column title
                "description": "The frequency of this barcode in the dataset",
                "min": 0,  # Maximum value in range, for bar / colour coding
                "format": "{:,." + decimal_positions + "f}",
            }

            columns_headers["Samples"] = {
                "namespace": "mgikit",
                "title": "Samples",  # Short title, table column title
                "description": "Samples that could match with this barcode",
            }
            self.add_section(
                name="Ambiguous barcodes by lane",
                anchor="mgikit-ambiguous-barcode",
                description="This Table shows the ambiguous barcodes that could match to multiple samples.",
                helptext="This longer string (can be **markdown**) helps explain how to interpret the plot",
                plot=table.plot(
                    sorted_dic,
                    columns_headers,
                    {
                        "id": "mgi_lane_ambiguous_plot",
                        "title": module_name + ": Ambiguous barcodes by lane",
                        "sort_rows": False,
                        "col1_header": "Barcode",
                    },
                ),
            )
        # else:
        #    log.info("No log files (*.mgikit.ambiguous_barcode) was found!")
