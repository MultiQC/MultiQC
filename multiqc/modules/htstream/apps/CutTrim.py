from collections import OrderedDict
import logging

from multiqc.plots import bargraph

#################################################

""" CutTrim submodule for HTStream charts and graphs """

#################################################


class CutTrim:
    ########################
    # Info about App
    def __init__(self):
        self.info = "Trims a fixed number of bases from the 5' and/or 3' end of each read."
        self.type = "bp_reducer"

    # Bargraph Function
    def bargraph(self, json, index):
        # config dict for bar graph
        config = {
            "title": "HTStream: Trimmed Basepairs Bargraph",
            "id": "htstream_cuttrim_bargraph_" + index,
            "ylab": "Percentage of Total Basepairs",
            "cpswitch": False,
            "data_labels": [
                {"name": "Percentage of Total", "ylab": "Percentage of Total Basepairs"},
                {"name": "Raw Counts", "ylab": "Basepairs"},
            ],
        }

        perc_data = {}
        read_data = {}

        # Construct data for multidataset bargraph
        for key in json:
            perc_data[key] = {
                "Perc_Left_Trim": json[key]["Ct_%_Left_Trimmed"],
                "Perc_Left_Trimm": json[key]["Ct_%_Right_Trimmed"],
            }
            read_data[key] = {"Left_Trim": json[key]["Ct_Left_Trimmed"], "Right_Trim": json[key]["Ct_Right_Trimmed"]}

        # Create categories for multidataset bargraph
        cats = [OrderedDict(), OrderedDict()]
        cats[0]["Perc_Left_Trim"] = {"name": "Left Trimmed"}
        cats[0]["Perc_Right_Trim"] = {"name": "Right Trimmed"}
        cats[1]["Left_Trim"] = {"name": "Left Trimmed"}
        cats[1]["Right_Trim"] = {"name": "Right Trimmed"}

        return bargraph.plot([perc_data, read_data], cats, config), ""

    ########################
    # MainFunction
    def execute(self, json, index):
        stats_json = OrderedDict()
        overview_dict = {}

        for key in json.keys():
            total_bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]
            left_bp_lost = (
                json[key]["Paired_end"]["Read1"]["leftTrim"]
                + json[key]["Paired_end"]["Read2"]["leftTrim"]
                + json[key]["Single_end"]["leftTrim"]
            )
            right_bp_lost = (
                json[key]["Paired_end"]["Read1"]["rightTrim"]
                + json[key]["Paired_end"]["Read2"]["rightTrim"]
                + json[key]["Single_end"]["rightTrim"]
            )

            try:
                perc_left_bp_lost = (left_bp_lost / json[key]["Fragment"]["basepairs_in"]) * 100
                perc_right_bp_lost = (right_bp_lost / json[key]["Fragment"]["basepairs_in"]) * 100
                fraction_bp_lost = total_bp_lost / json[key]["Fragment"]["basepairs_in"]

            except ZeroDivisionError:
                perc_left_bp_lost = 0
                perc_right_bp_lost = 0
                fraction_bp_lost = 0

                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)
                continue

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Bp_Lost": fraction_bp_lost,
            }

            # sample dictionary entry
            stats_json[key] = {
                "Ct_%_Left_Trimmed": perc_left_bp_lost,
                "Ct_%_Right_Trimmed": perc_right_bp_lost,
                "Ct_Left_Trimmed": left_bp_lost,
                "Ct_Right_Trimmed": right_bp_lost,
            }

        # sections and figure function calls
        figure, html = self.bargraph(stats_json, index)
        section = {"Figure": figure, "Overview": overview_dict, "Content": html}

        return section
