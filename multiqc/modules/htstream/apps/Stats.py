from collections import OrderedDict
import logging
from random import random

from . import htstream_utils
from multiqc.plots import table, linegraph

#################################################

""" Stats submodule for HTStream charts and graphs """

#################################################


class Stats:
    ########################
    # Info about App
    def __init__(self):
        self.info = "Generates a summary report from a set of statistical measures about the input read data."
        self.type = "both"
        self.read_keys = {
            "St_PE_Base_by_Cycle": "PE",
            "St_SE_Base_by_Cycle": "SE",
            "St_PE_Quality_by_Cycle": "PE",
            "St_SE_Quality_by_Cycle": "SE",
            "St_PE_Read_Lengths": "PE",
            "St_SE_Read_Lengths": "SE",
            "PE": "Paired End",
            "SE": "Single End",
        }

    ########################
    # Table Function
    def table(self, json, index):
        unique_id = str(random())[2:6]

        # pconfig
        config = {"id": "htstream_stats_table_line_" + unique_id, "title": "HTStream: Stats"}

        # striaght forward table function, right from MultiQC documentation
        headers = OrderedDict()

        # "St_PE_Fraction" + index
        headers["St_R1_Q30" + index] = {
            "title": "% R1 Q30",
            "namespace": "% R1 Q30",
            "description": "percentage of read 1 bps Q30 or greater",
            "format": "{:,.2f}",
            "suffix": "%",
            "scale": "RdBu",
        }
        headers["St_R2_Q30" + index] = {
            "title": "% R2 Q30",
            "namespace": "% R2 Q30",
            "description": "percentage of read 2 bps Q30 or greater",
            "format": "{:,.2f}",
            "suffix": "%",
            "scale": "PRGn",
        }
        headers["St_SE_Q30" + index] = {
            "title": "% SE Q30",
            "namespace": "% SE Q30",
            "description": "percentage of single end read bps Q30 or greater",
            "format": "{:,.2f}",
            "suffix": "%",
            "scale": "RdPu",
        }
        headers["St_GC_Content" + index] = {
            "title": "GC Content",
            "namespace": "GC Content",
            "description": "Percentage of bps that are G or C",
            "format": "{:,.2f}",
            "suffix": "%",
            "scale": "Reds",
        }
        headers["St_N_Content" + index] = {
            "title": "N Content",
            "namespace": "N Content",
            "description": "Percentage of bps that are N",
            "format": "{:,.4f}",
            "suffix": "%",
            "scale": "RdPu",
        }

        return table.plot(json, headers, config)

    ########################
    # Base By cycle sections Functions
    def base_by_cycle(self, json, read):
        # Read Code and Unique ID
        read_code = self.read_keys[read]
        unique_id = str(random())[2:6]

        # Multi Sample Line Config, it's called entropy (even though its not)
        config = {
            "id": "htstream_stats_distance_" + read + "_" + unique_id,
            "title": "HTStream: Base by Cycle (" + read_code + ")",
            "smooth_points_sumcounts": False,
            "ylab": "Avg. Distance from 25%",
            "xlab": "Cycle",
            "categories": True,
            "tt_label": "{point.x}: {point.y:.2f}%",
            "y_bands": [
                {"from": 0, "to": 8, "color": "#c3e6c3"},
                {"from": 8, "to": 35, "color": "#e6dcc3"},
                {"from": 35, "to": 100, "color": "#e6c3c3"},
            ],
        }

        # If paired end, we need to add midpoint line
        if read_code == "PE":
            midpoint = htstream_utils.uniform(json, read)

            # If uniform, add line
            if midpoint != -1:
                line_list = [
                    {
                        "color": "#5D4B87",
                        "width": 1,
                        "value": (midpoint - 1) / 2,
                        "dashStyle": "shortdash",
                    }
                ]

                config["x_lines"] = line_list

        # Data list and dictionaries for line graphs
        line_data = {}

        # For each sample, add their line to multi-sample graph and construct their own line graph
        for samp in json.keys():
            line_data[samp] = {}

            # If PE< we need to concat base by cycle lists
            if read_code == "PE":
                data = [json[samp][read][0]["data"][x] + json[samp][read][1]["data"][x] for x in range(5)]
            else:
                data = json[samp][read]["data"]

            # Len of reads
            bases = len(data[0])

            # Iterates through positions
            for i in range(bases):
                # Total count at position
                total = sum([base[i] for base in data])

                # Avg difference from 25%, N not included
                avg = sum([abs(((data[x][i] / total) * 100) - 25) for x in range(4)]) / 4

                line_data[samp][i + 1] = avg

        return linegraph.plot(line_data, config)

    ########################
    # Quality By cycle sections Functions
    def quality_by_cycle(self, json, read):
        # Read Code and Unique ID
        read_code = self.read_keys[read]
        unique_id = str(random())[2:6]

        # config dictionary for mean Q score line graph
        line_config = {
            "id": "htstream_stats_qbc_line_" + read_code + "_" + unique_id,
            "smooth_points_sumcounts": False,
            "categories": True,
            "tt_label": "{point.x}: {point.y:.2f}",
            "title": "HTStream: Mean Quality by Cycle (" + read_code + ")",
            "xlab": "Cycle",
            "ylab": "Mean Q Score",
            "colors": {},
        }

        # If paired end, we need to add midpoint line
        if read_code == "PE":
            midpoint = htstream_utils.uniform(json, read)

            # If uniform, add line
            if midpoint != -1:
                line_config["x_lines"] = [
                    {
                        "color": "#5D4B87",
                        "width": 1,
                        "value": (midpoint - 1) / 2,
                        "dashStyle": "shortdash",
                    }
                ]

        # Line Data, Button List, and Boolean
        line_data = {}

        for key in json.keys():
            # create dictionary for line graph. Again, format is {x: y}
            line_data[key] = {}

            # If PE, we wanna concat quality by cycle data
            if read_code == "PE":
                # Length of R1, Create concat list, Get Column Names
                length = len(json[key][read][0]["data"])
                temp_data = [json[key][read][0]["data"][x] + json[key][read][1]["data"][x] for x in range(length)]
                temp_col_name = [
                    str(int(x) + int(json[key][read][0]["col_names"][-1])) for x in json[key][read][1]["col_names"]
                ]

                # Rewrite JSON for this sample
                json[key][read] = {
                    "data": temp_data,
                    "col_names": json[key][read][0]["col_names"] + temp_col_name,
                    "row_names": json[key][read][0]["row_names"],
                    "shape": [
                        json[key][read][0]["shape"][0],
                        json[key][read][0]["shape"][-1] + json[key][read][1]["shape"][-1],
                    ],
                }

            y_lab = json[key][read]["row_names"][::-1]  # reverse orientation makes it easier to cycle through

            # create variables for range functions in loops. Represents shape of data
            # quality_scores = json[key][read]["shape"][0]
            cycles = json[key][read]["shape"][-1]

            # iterates through positions, creates a list of the sum of scores at each position to be used
            # 	to calculated frequency for heatmap. Also, calculates avg. Q score for linegraph.
            # 	This chunk of code is very ugly, but is a necessary evil.

            num_above_q30 = 0

            for pos in range(cycles):
                temp = [score_list[pos] for score_list in json[key][read]["data"]]
                temp_sum = sum(temp)

                # multiples count at poistion by Q Score.
                total_score = sum([(int(p) * int(s)) for p, s in zip(temp, y_lab[::-1])])

                # divides sum of total score by the number of cycles for avg fragments
                line_data[key][pos + 1] = total_score / temp_sum  # total reads

                if line_data[key][pos + 1] > 30:
                    num_above_q30 += 1

            # check to see what percent of bases have a mean Q score of at least 30,
            #   color samples accordingly.
            q30_gate = num_above_q30 / cycles

            if q30_gate < 0.6:
                line_config["colors"][key] = "#E16B6B"

            elif q30_gate < 0.8:
                line_config["colors"][key] = "#E8A243"

            else:
                line_config["colors"][key] = "#78D578"

        return linegraph.plot(line_data, line_config)

    ########################
    # Read Length Heatmaps
    def read_length(self, json, read):
        # Read cor and unique IDs
        read_code = self.read_keys[read]
        unique_id = str(random())[2:6]

        # config dictionary for linegraph
        line_config = {
            "id": "htstream_stats_read_lengths_" + read_code + "_" + unique_id,
            "title": "HTStream: Read Length Linegraph (" + read_code + ")",
            "ylab": "Counts",
            "xlab": "Read Length",
        }

        # if PE, cread multiple datasets
        if read_code == "PE":
            line_config["data_labels"] = [
                {"name": "Read 1", "ylab": "Counts", "xlab": "Read Lengths"},
                {"name": "Read 2", "ylab": "Counts", "xlab": "Read Lengths"},
            ]
            readlength_data = [{}, {}]

        else:
            readlength_data = {}

        # Are all reads from all samples uniform? Let's find out.
        uniform_dict = {}

        for samp in json.keys():
            # paired end reads require the histograms be concatenated
            if read_code == "SE":
                # uniform reads
                if len(json[samp][read][0]) == 1:
                    uniform_dict[samp] = {"St_Read_Lengths_SE_" + unique_id: json[samp][read][0][0][0]}

                # create entry
                readlength_data[samp] = {}

                # populate x,y coords
                for length in json[samp][read][0]:
                    readlength_data[samp][length[0]] = length[1]

            else:
                # uniform read lengths
                if len(json[samp][read][0]) + len(json[samp][read][1]) == 2:
                    uniform_dict[samp] = {
                        "St_Read_Lengths_R1_" + unique_id: json[samp][read][0][0][0],
                        "St_Read_Lengths_R2_" + unique_id: json[samp][read][1][0][0],
                    }

                # iterate through R1 and R2 read length hists
                for i in range(0, 2):
                    readlength_data[i][samp] = {}

                    for length in json[samp][read][i]:
                        readlength_data[i][samp][length[0]] = length[1]

        if len(uniform_dict.keys()) == len(json.keys()):
            config = {"id": "htstream_stats_read_length_" + unique_id, "title": "HTStream: Stats"}

            headers = OrderedDict()

            headers["St_Read_Lengths_SE_" + unique_id] = {
                "title": "SE Read Lengths",
                "namespace": "SE Read Lengths",
                "description": "Length of SE reads (uniform for each sample).",
                "format": "{:,.0f}",
                "scale": "Blues",
            }
            headers["St_Read_Lengths_R1_" + unique_id] = {
                "title": "R1 Read Lengths",
                "namespace": "R1 Read Lengths",
                "description": "Length of R1 reads (uniform for each sample).",
                "format": "{:,.0f}",
                "scale": "Blues",
            }
            headers["St_Read_Lengths_R2_" + unique_id] = {
                "title": "R2 Read Lengths",
                "namespace": "R2 Read Lengths",
                "description": "Length of R2 reads (uniform for each sample).",
                "format": "{:,.0f}",
                "scale": "Greens",
            }

            figure = table.plot(uniform_dict, headers, config)

        else:
            figure = linegraph.plot(readlength_data, line_config)

        return figure

    ########################
    # Main Function
    def execute(self, json, index):
        stats_json = OrderedDict()
        SE_json = OrderedDict()
        PE_json = OrderedDict()
        overview_stats = {}

        for key in json.keys():
            #
            # STATS FOR TABLE
            #
            total_frags = json[key]["Fragment"]["out"]
            total_bps = json[key]["Fragment"]["basepairs_out"]

            # exit if no reads, report that something is wrong
            if total_frags == 0 or total_bps == 0:
                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)
                continue

            gc_count = json[key]["Fragment"]["base_composition"]["G"] + json[key]["Fragment"]["base_composition"]["C"]
            gc_content = gc_count / total_bps
            n_content = json[key]["Fragment"]["base_composition"]["N"] / total_bps

            stats_json[key] = {
                "St_Fragments_in" + index: total_frags,
                "St_GC_Content" + index: gc_content * 100,
                "St_N_Content" + index: n_content * 100,
            }

            # Overview stats
            overview_stats[key] = {
                "GC_Content": gc_content,
                "N_Content": n_content,
                "Q30_Fraction": 0,
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
            }
            #
            # SINGLE END STATS
            #
            # only succeeds if json file contains single end information data in the last instance of hts_Stats,
            # 	opens gate for future processing of single end read stats.
            try:
                SE_json[key] = {}
                SE_json[key]["St_SE_Read_Lengths"] = [json[key]["Single_end"]["readlength_histogram"]]
                SE_json[key]["St_SE_Base_by_Cycle"] = json[key]["Single_end"]["base_by_cycle"]
                SE_json[key]["St_SE_Quality_by_Cycle"] = json[key]["Single_end"]["qualities_by_cycle"]
                SE_json[key]["St_SE_in"] = json[key]["Single_end"]["in"]

                stats_json[key]["St_SE_Q30" + index] = (
                    json[key]["Single_end"]["total_Q30_basepairs"] / json[key]["Single_end"]["basepairs_in"]
                ) * 100

                overview_stats[key]["Q30_Fraction"] += json[key]["Single_end"]["total_Q30_basepairs"] / total_bps
                overview_stats[key]["SE_Fraction"] = json[key]["Single_end"]["out"] / total_frags

            except KeyError:
                overview_stats[key]["Q30_Fraction"] += 0
                overview_stats[key]["SE_Fraction"] = 0

                del SE_json[key]

            #
            # PAIRED END STATS
            #
            try:
                PE_json[key] = {}
                PE_json[key]["St_PE_Read_Lengths"] = [
                    json[key]["Paired_end"]["Read1"]["readlength_histogram"],
                    json[key]["Paired_end"]["Read2"]["readlength_histogram"],
                ]

                PE_json[key]["St_PE_Base_by_Cycle"] = [
                    json[key]["Paired_end"]["Read1"]["base_by_cycle"],
                    json[key]["Paired_end"]["Read2"]["base_by_cycle"],
                ]
                PE_json[key]["St_PE_Quality_by_Cycle"] = [
                    json[key]["Paired_end"]["Read1"]["qualities_by_cycle"],
                    json[key]["Paired_end"]["Read2"]["qualities_by_cycle"],
                ]
                PE_json[key]["St_PE_in"] = json[key]["Paired_end"]["in"]

                stats_json[key]["St_R1_Q30" + index] = (
                    json[key]["Paired_end"]["Read1"]["total_Q30_basepairs"]
                    / json[key]["Paired_end"]["Read1"]["basepairs_in"]
                ) * 100
                stats_json[key]["St_R2_Q30" + index] = (
                    json[key]["Paired_end"]["Read2"]["total_Q30_basepairs"]
                    / json[key]["Paired_end"]["Read2"]["basepairs_in"]
                ) * 100

                overview_stats[key]["Q30_Fraction"] += (
                    json[key]["Paired_end"]["Read1"]["total_Q30_basepairs"]
                    + json[key]["Paired_end"]["Read2"]["total_Q30_basepairs"]
                ) / total_bps
                overview_stats[key]["PE_Fraction"] = json[key]["Paired_end"]["out"] / total_frags

            except KeyError:
                overview_stats[key]["Q30_Fraction"] += 0
                overview_stats[key]["PE_Fraction"] = 0

                del PE_json[key]

        # output dictionary, keys are section, value is function called for figure generation
        section = {"Table": self.table(stats_json, index), "Overview": overview_stats}

        if len(PE_json.keys()) != 0:
            section["PE"] = {
                "Read_Length": self.read_length(PE_json, "St_PE_Read_Lengths"),
                "Base_by_Cycle": self.base_by_cycle(PE_json, "St_PE_Base_by_Cycle"),
                "Quality_by_Cycle": self.quality_by_cycle(PE_json, "St_PE_Quality_by_Cycle"),
            }

        # only executres if single read data is detected
        if len(SE_json.keys()) != 0:
            section["SE"] = {
                "Read_Length": self.read_length(SE_json, "St_SE_Read_Lengths"),
                "Base_by_Cycle": self.base_by_cycle(SE_json, "St_SE_Base_by_Cycle"),
                "Quality_by_Cycle": self.quality_by_cycle(SE_json, "St_SE_Quality_by_Cycle"),
            }

        return section
