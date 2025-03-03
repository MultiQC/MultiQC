from collections import OrderedDict
import logging

from multiqc.plots import bargraph

#################################################

""" LengthFilter submodule for HTStream charts and graphs """

#################################################


class LengthFilter:
    ########################
    # Info about App
    def __init__(self):
        self.info = "Discards reads below a minimum length threshold."
        self.type = "read_reducer"

    # Bargraph Function
    def bargraph(self, json, reads_lost, index):
        # config dict for bar graph
        config = {
            "title": "HTStream: Removed Reads Bargraph",
            "id": "htstream_lengthfilter_bargraph_" + index,
            "cpswitch": False,
            "ylab": "Percentage of Total Reads",
            "data_labels": [
                {"name": "Percentage of Total", "ylab": "Percentage of Total Reads"},
                {"name": "Raw Counts", "ylab": "Reads"},
            ],
        }

        html = ""

        # returns nothing if no reads were trimmed.
        if reads_lost == 0:
            html += (
                '<div class="alert alert-info"> <strong>Notice:</strong> No reads were removed in any sample. </div>'
            )
            figure = None

        else:
            perc_data = {}
            read_data = {}

            # Construct data for multidataset bargraph
            for key in json:
                perc_pe = (json[key]["Lf_PE_lost"] / json[key]["Lf_Total_Reads"]) * 100
                perc_se = (json[key]["Lf_SE_lost"] / json[key]["Lf_Total_Reads"]) * 100

                perc_data[key] = {"Perc_PE": perc_pe, "Perc_SE": perc_se}
                read_data[key] = {
                    "Reads_PE": json[key]["Lf_PE_lost"],
                    "Reads_SE": json[key]["Lf_SE_lost"],
                }

            # Create categories for multidataset bargraph
            cats = [OrderedDict(), OrderedDict()]
            cats[0]["Perc_PE"] = {"name": "Paired End"}
            cats[0]["Perc_SE"] = {"name": "Single End"}
            cats[1]["Reads_PE"] = {"name": "Paired End"}
            cats[1]["Reads_SE"] = {"name": "Single End"}

            # Create bargraph
            figure = bargraph.plot([perc_data, read_data], cats, config)

        return figure, html

    ########################
    # Main Function
    def execute(self, json, index):
        stats_json = OrderedDict()
        overview_dict = {}

        # Accumulator vars
        reads_lost = 0

        for key in json.keys():
            if json[key]["Fragment"]["in"] == 0:
                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)

            reads_lost += json[key]["Fragment"]["in"] - json[key]["Fragment"]["out"]

            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Reads_Lost": reads_lost,
            }

            # sample entry for stats dictionary
            stats_json[key] = {
                "Lf_Total_Reads": json[key]["Fragment"]["in"],
                "Lf_PE_lost": json[key]["Paired_end"]["discarded"],
                "Lf_SE_lost": json[key]["Single_end"]["discarded"],
            }

        figure, html = self.bargraph(stats_json, reads_lost, index)

        # sections and figure function calls
        section = {"Figure": figure, "Overview": overview_dict, "Content": html}

        return section
