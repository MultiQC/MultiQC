from collections import OrderedDict
import logging

from multiqc.plots import bargraph

#################################################

""" QWindowTrim submodule for HTStream charts and graphs """

#################################################


class QWindowTrim:
    ########################
    # Info about App
    def __init__(self):
        self.info = "Uses a sliding window approach to remove the low quality ends of reads."
        self.type = "bp_reducer"

    ########################
    # Table Function
    def bargraph(self, json, bps, index):
        # config dict for bar graph
        config = {
            "title": "HTStream: QWindowTrim Trimmed Basepairs Bargraph",
            "id": "htstream_qwindowtrim_bargraph_" + index,
            "ylab": "Percentage of Total Basepairs",
            "cpswitch": False,
            "data_labels": [
                {"name": "Read 1", "ylab": "Percentage of Total Basepairs"},
                {"name": "Read 2", "ylab": "Percentage of Total Basepairs"},
                {"name": "Single End", "ylab": "Percentage of Total Basepairs"},
            ],
        }

        # Header
        html = ""

        # returns nothing if no reads were trimmed.
        if bps == 0:
            html += '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from any sample. </div>'
            figure = None

        else:
            r1_data = {}
            r2_data = {}
            se_data = {}

            # Create dictionaries for multidataset bargraphs
            for key in json:
                r1_data[key] = {"LT_R1": json[key]["Qt_Left_Trimmed_R1"], "RT_R1": json[key]["Qt_Right_Trimmed_R1"]}

                r2_data[key] = {"LT_R2": json[key]["Qt_Left_Trimmed_R2"], "RT_R2": json[key]["Qt_Right_Trimmed_R2"]}

                se_data[key] = {"LT_SE": json[key]["Qt_Left_Trimmed_SE"], "RT_SE": json[key]["Qt_Right_Trimmed_SE"]}

            # Create categores for multidatatset bragraphs
            cats = [OrderedDict(), OrderedDict(), OrderedDict()]
            cats[0]["LT_R1"] = {"name": "Left Trimmmed"}
            cats[0]["RT_R1"] = {"name": "Right Trimmmed"}
            cats[1]["LT_R2"] = {"name": "Left Trimmmed"}
            cats[1]["RT_R2"] = {"name": "Right Trimmmed"}
            cats[2]["LT_SE"] = {"name": "Left Trimmmed"}
            cats[2]["RT_SE"] = {"name": "Right Trimmmed"}

            figure = bargraph.plot([r1_data, r2_data, se_data], cats, config)

        return figure, html

    ########################
    # Main Function
    def execute(self, json, index):
        stats_json = OrderedDict()
        overview_dict = {}

        overall_trim = 0

        for key in json.keys():
            total_bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]
            overall_trim += total_bp_lost

            bp_in = json[key]["Fragment"]["basepairs_in"]

            if bp_in != 0:
                fract_r1_bp_left = json[key]["Paired_end"]["Read1"]["leftTrim"] / bp_in
                fract_r1_bp_right = json[key]["Paired_end"]["Read1"]["rightTrim"] / bp_in
                fract_r2_bp_left = json[key]["Paired_end"]["Read2"]["leftTrim"] / bp_in
                fract_r2_bp_right = json[key]["Paired_end"]["Read2"]["rightTrim"] / bp_in
                fract_se_bp_left = json[key]["Single_end"]["leftTrim"] / bp_in
                fract_se_bp_right = json[key]["Single_end"]["rightTrim"] / bp_in

            else:
                fract_r1_bp_left = 0
                fract_r1_bp_right = 0
                fract_r2_bp_left = 0
                fract_r2_bp_right = 0
                fract_se_bp_left = 0
                fract_se_bp_right = 0

                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)

            # overview data
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_R1_Bp_Trimmed_Left": fract_r1_bp_left,
                "Fraction_R1_Bp_Trimmed_Right": fract_r1_bp_right,
                "Fraction_R2_Bp_Trimmed_Left": fract_r2_bp_left,
                "Fraction_R2_Bp_Trimmed_Right": fract_r2_bp_right,
                "Fraction_SE_Bp_Trimmed_Left": fract_se_bp_left,
                "Fraction_SE_Bp_Trimmed_Right": fract_se_bp_right,
            }

            # sample dictionary entry
            stats_json[key] = {
                "Qt_Left_Trimmed_R1": fract_r1_bp_left * 100,
                "Qt_Right_Trimmed_R1": fract_r1_bp_right * 100,
                "Qt_Left_Trimmed_R2": fract_r2_bp_left * 100,
                "Qt_Right_Trimmed_R2": fract_r2_bp_right * 100,
                "Qt_Left_Trimmed_SE": fract_se_bp_left * 100,
                "Qt_Right_Trimmed_SE": fract_se_bp_right * 100,
            }

        # sections and figure function calls
        figure, html = self.bargraph(stats_json, overall_trim, index)
        section = {"Figure": figure, "Overview": overview_dict, "Content": html}

        return section
