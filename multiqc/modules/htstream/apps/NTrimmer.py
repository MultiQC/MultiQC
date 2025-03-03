from collections import OrderedDict
import logging

from multiqc.plots import bargraph

#################################################

""" NTrimmer submodule for HTStream charts and graphs """

#################################################


class NTrimmer:
    ########################
    # Info about App
    def __init__(self):
        self.info = "Trims reads to the longest subsequence that contains no N's."
        self.type = "bp_reducer"

    ########################
    # Table Function
    def bargraph(self, json, bps, index):
        # config dict for bar graph
        config = {
            "title": "HTStream: NTrimmer Trimmed Basepairs Bargraph",
            "id": "htstream_ntrimmer_bargraph_" + index,
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
            html = '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from any sample. </div>'
            figure = None

        else:
            r1_data = {}
            r2_data = {}
            se_data = {}

            # Create dictionaries for multidataset bargraphs
            for key in json:
                r1_data[key] = {"LT_R1": json[key]["Nt_Left_Trimmed_R1"], "RT_R1": json[key]["Nt_Right_Trimmed_R1"]}

                r2_data[key] = {"LT_R2": json[key]["Nt_Left_Trimmed_R2"], "RT_R2": json[key]["Nt_Right_Trimmed_R2"]}

                se_data[key] = {"LT_SE": json[key]["Nt_Left_Trimmed_SE"], "RT_SE": json[key]["Nt_Right_Trimmed_SE"]}

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

        # accumulator variable. Used to prevent empty bargraphs
        overall = 0

        for key in json.keys():
            total_bp = json[key]["Fragment"]["basepairs_in"]
            total_bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]

            total_r1 = (
                json[key]["Paired_end"]["Read1"]["basepairs_in"] - json[key]["Paired_end"]["Read1"]["basepairs_out"]
            )

            total_r2 = (
                json[key]["Paired_end"]["Read2"]["basepairs_in"] - json[key]["Paired_end"]["Read2"]["basepairs_out"]
            )
            total_se = json[key]["Single_end"]["basepairs_in"] - json[key]["Single_end"]["basepairs_out"]

            try:
                fract_bp_lost = total_bp_lost / json[key]["Fragment"]["basepairs_in"]

            except ZeroDivisionError:
                fract_bp_lost = 0

                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)
                raise

            # overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Bp_Lost": fract_bp_lost,
            }

            # sample entry in stats dictionary
            stats_json[key] = {
                "Nt_Left_Trimmed_R1": (json[key]["Paired_end"]["Read1"]["leftTrim"] / total_bp) * 100,
                "Nt_Right_Trimmed_R1": (json[key]["Paired_end"]["Read1"]["rightTrim"] / total_bp) * 100,
                "Nt_Left_Trimmed_R2": (json[key]["Paired_end"]["Read2"]["leftTrim"] / total_bp) * 100,
                "Nt_Right_Trimmed_R2": (json[key]["Paired_end"]["Read2"]["rightTrim"] / total_bp) * 100,
                "Nt_Left_Trimmed_SE": (json[key]["Single_end"]["leftTrim"] / total_bp) * 100,
                "Nt_Right_Trimmed_SE": (json[key]["Single_end"]["rightTrim"] / total_bp) * 100,
            }

            # accumulatr totals
            overall += total_r1 + total_r2 + total_se

        # section and figure function calls
        figure, html = self.bargraph(stats_json, overall, index)
        section = {"Figure": figure, "Overview": overview_dict, "Content": html}

        return section
