from collections import OrderedDict
import logging

from multiqc.plots import bargraph

#################################################

""" Primers submodule for HTStream charts and graphs """

#################################################


class Primers:
    ########################
    # Info about App
    def __init__(self):
        self.info = "Identifies primer sequences located on the 5' ends of R1 and R2, or 5' and 3' end of SE reads."
        self.type = "bp_reducer"

    # Bargraph Function
    def bargraph(self, json, index):
        # config dict for bar graph
        config = {
            "title": "HTStream: Primer Counts Bargraph",
            "id": "htstream_primers_bargraph_" + index,
            "name": "Primer Counts",
            "ylab": "Primer Combination Counts",
        }

        html = ""

        data = {}
        labels = []

        # Construct data for multidataset bargraph
        for key in json:
            counts_list = json[key]["Pr_Primer_Counts"]
            data[key] = {}

            # get counts and labels
            for x in range(len(counts_list)):
                label = "-".join(counts_list[x][:-1])
                data[key][label] = counts_list[x][-1]
                labels.append(label)

        cats = list(set(labels))

        return bargraph.plot(data, cats, config), html

    ########################
    # Main Function
    def execute(self, json, index):
        stats_json = OrderedDict()
        overview_dict = {}

        for key in json.keys():
            bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]

            try:
                fract_bp_lost = bp_lost / json[key]["Fragment"]["basepairs_in"]

            except ZeroDivisionError:
                fract_bp_lost = 0

                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)
                raise

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Bp_Lost": fract_bp_lost,
            }

            stats_json[key] = {"Pr_Primer_Counts": json[key]["Fragment"]["primers_counts"]}

        # dictionary for sections and figure function calls
        figure, html = self.bargraph(stats_json, index)
        section = {"Figure": figure, "Overview": overview_dict, "Content": html}

        return section
