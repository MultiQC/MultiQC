from collections import OrderedDict
import logging

from multiqc.plots import bargraph

#################################################

""" PolyATTrim submodule for HTStream charts and graphs """

#################################################


class PolyATTrim:
    ########################
    # Info about App
    def __init__(self):
        self.info = "Attempts to trim poly-A and poly-T sequences from the end of reads."
        self.type = "bp_reducer"

    ########################
    # Bargraph Function
    def bargraph(self, json, bps_trimmed, index):
        # configuration dictionary for bar graph
        config = {
            "title": "HTStream: PolyATTrim Trimmed Basepairs Bargraph",
            "id": "htstream_polyattrim_bargraph_" + index,
            "ylab": "Percentage of Total Basepairs",
            "cpswitch": False,
            "data_labels": [
                {"name": "Percentage of Total", "ylab": "Percentage of Total Basepairs"},
                {"name": "Raw Counts", "ylab": "Basepairs"},
            ],
        }

        # Title
        html = ""

        # if no overlaps at all are present, return nothing
        if bps_trimmed == 0:
            html += (
                '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from samples. </div>'
            )
            figure = None

        else:
            perc_data = {}
            read_data = {}

            # Construct data for multidataset bargraph
            for key in json:
                perc_data[key] = {
                    "Perc_R1_lost": json[key]["Pt_Perc_R1_lost"],
                    "Perc_R2_lost": json[key]["Pt_Perc_R2_lost"],
                    "Perc_SE_lost": json[key]["Pt_Perc_SE_lost"],
                }
                read_data[key] = {
                    "R1_lost": json[key]["Pt_R1_lost"],
                    "R2_lost": json[key]["Pt_R2_lost"],
                    "SE_lost": json[key]["Pt_SE_lost"],
                }

            # bargraph dictionary. Exact use of example in MultiQC docs.
            categories = [OrderedDict(), OrderedDict()]

            # Colors for sections
            categories[0]["Perc_R1_lost"] = {"name": "Read 1", "color": "#779BCC"}
            categories[0]["Perc_R2_lost"] = {"name": "Read 2", "color": "#C3C3C3"}
            categories[0]["Perc_SE_lost"] = {"name": "Single End", "color": "#D1ADC3"}
            categories[1]["R1_lost"] = {"name": "Read 1", "color": "#779BCC"}
            categories[1]["R2_lost"] = {"name": "Read 2", "color": "#C3C3C3"}
            categories[1]["SE_lost"] = {"name": "Single End", "color": "#D1ADC3"}

            figure = bargraph.plot([perc_data, read_data], categories, config)

        return figure, html

    ########################
    # Main Function
    def execute(self, json, index):
        stats_json = OrderedDict()
        overview_dict = {}

        # accumulator variable. Used to prevent empty bargraphs
        overall = 0

        for key in json.keys():
            total_bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]

            r1_lost = (
                json[key]["Paired_end"]["Read1"]["basepairs_in"] - json[key]["Paired_end"]["Read1"]["basepairs_out"]
            )
            r2_lost = (
                json[key]["Paired_end"]["Read2"]["basepairs_in"] - json[key]["Paired_end"]["Read2"]["basepairs_out"]
            )
            se_lost = json[key]["Single_end"]["basepairs_in"] - json[key]["Single_end"]["basepairs_out"]

            # try block to avoid zero division
            try:
                fract_bp_lost = total_bp_lost / json[key]["Fragment"]["basepairs_in"]

                perc_r1_lost = (r1_lost / json[key]["Fragment"]["basepairs_in"]) * 100
                perc_r2_lost = (r2_lost / json[key]["Fragment"]["basepairs_in"]) * 100
                perc_se_lost = (se_lost / json[key]["Fragment"]["basepairs_in"]) * 100

            except ZeroDivisionError:
                fract_bp_lost = 0

                perc_r1_lost = 0
                perc_r2_lost = 0
                perc_se_lost = 0

                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)
                continue

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Bp_Lost": fract_bp_lost,
            }

            # sample entry in stats dictionary
            stats_json[key] = {
                "Pt_Perc_R1_lost": perc_r1_lost,
                "Pt_Perc_R2_lost": perc_r2_lost,
                "Pt_Perc_SE_lost": perc_se_lost,
                "Pt_R1_lost": r1_lost,
                "Pt_R2_lost": r2_lost,
                "Pt_SE_lost": se_lost,
            }

            # Accumulate totals
            overall += total_bp_lost

        # section and figure function calls
        figure, html = self.bargraph(stats_json, overall, index)
        section = {"Figure": figure, "Overview": overview_dict, "Content": html}

        return section
