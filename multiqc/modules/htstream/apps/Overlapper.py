from collections import OrderedDict
import logging
import statistics

from multiqc.plots import bargraph

#################################################

""" Overlapper submodule for HTStream charts and graphs """

#################################################


class Overlapper:
    ########################
    # Info about App
    def __init__(self):
        self.info = "Attempts to overlap paired end reads to produce the original fragment, trims adapters, and can correct sequencing errors."
        self.type = "read_reducer"

    ########################
    # Bargraph Function
    def bargraph(self, json, inserts):
        # configuration dictionary for bar graph
        config = {
            "title": "HTStream: Overlap Composition Bargraph",
            "id": "htstream_overlapper_bargraph",
            "ylab": "Overlap Types",
            "cpswitch": True,
            "cpswitch_c_active": True,
            "cpswitch_counts_label": "Counts",
            "cpswitch_percent_label": "Percentages",
        }

        # Header
        html = ""

        # if no overlaps at all are present, return nothing
        if inserts == 0:
            html = '<div class="alert alert-info"> <strong>Notice:</strong> No overlaps present in samples. </div>'
            figure = None

        else:
            # bargraph dictionary. Exact use of example in MultiQC docs.
            categories = OrderedDict()

            # Create blocks for bargrapph
            categories["Ov_Sins"] = {"name": "Short Inserts", "color": "#779BCC"}
            categories["Ov_Mins"] = {"name": "Medium Inserts", "color": "#C3C3C3"}
            categories["Ov_Lins"] = {"name": "Long Inserts", "color": "#D1ADC3"}

            # create plot
            figure = bargraph.plot(json, categories, config)

        return figure, html

    # ########################
    # # Linegraph Function
    # def linegraph(self, json, index):

    #     # config dictionary for "density" plots. Its a work in progress.
    #     config = {
    #         "id": "htstream_overlapper_linegraph_" + index,
    #         "title": "HTStream: Overlapped Lengths",
    #         "ylab": "Counts",
    #         "xlab": "Overlap Lengths",
    #     }

    #     # initialize data structures
    #     multi_line = {}

    #     for key in json.keys():

    #         # creates empty dictionary to hold data for line graph.
    #         multi_line[key] = {}

    #         # iterates over ever value in histogram and adds it to line graph
    #         for item in json[key]["Ov_Histogram"]:

    #             multi_line[key][item[0]] = item[1]

    #     html = "<h4> Overlapper: Overlapped Lengths </h4>\n"
    #     html += "<p>Plots the lengths of paired end read overlaps.</p>"
    #     html += linegraph.plot(multi_line, config)

    #     return html

    ########################
    # Function for parsing histogram in files
    def parse_histogram_stats(self, hist):
        hist_stats = {"Max": 0, "Median": 0}

        median_list = []

        # Find median and max of list in files
        for item in hist:
            median_list.append(item[0])

            if hist_stats["Max"] < item[1]:
                hist_stats["Max"] = item[0]

        hist_stats["Median"] = statistics.median(median_list)

        return hist_stats

    ########################
    # Main Function
    def execute(self, json, index):
        stats_json = OrderedDict()
        overview_dict = {}

        # accumulator for inserts, used to prevent empty bar graph
        inserts = 0
        # se_total_gain = 0

        for key in json.keys():
            sins = json[key]["Fragment"]["inserts"]["short"]
            mins = json[key]["Fragment"]["inserts"]["medium"]
            lins = json[key]["Fragment"]["inserts"]["long"]

            # Total overlap types
            overlapped_sum = sins + mins + lins

            # try to parse overlap hist
            try:
                ov_hist = json[key]["Fragment"].get("overlap_histogram", [[0, 0]])
                parsed_hist_stats = self.parse_histogram_stats(ov_hist)
            except:
                parsed_hist_stats = -1
                ov_hist = -1
                raise

            # if no histogram, assign zeroes to median and max
            if parsed_hist_stats == -1:
                ov_max = 0
                ov_med = 0

            else:
                ov_max = parsed_hist_stats["Max"]
                ov_med = parsed_hist_stats["Median"]

            # try for perc of sins, mins, and lins,
            try:
                perc_sin = sins / json[key]["Fragment"]["in"]
                perc_min = mins / json[key]["Fragment"]["in"]
                perc_lin = lins / json[key]["Fragment"]["in"]

            except:
                perc_sin = 0
                perc_min = 0
                perc_lin = 0

                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)
                raise

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Overlap_Length_Max": ov_max,
                "Overlap_Length_Med": ov_med,
                "Sin": perc_sin,
                "Min": perc_min,
                "Lin": perc_lin,
            }

            # sample instance in dictionary
            stats_json[key] = {"Ov_Sins": sins, "Ov_Mins": mins, "Ov_Lins": lins, "Ov_Histogram": ov_hist}

            # accumulator accumlating
            inserts += overlapped_sum

        # sections and function calls
        figure, html = self.bargraph(stats_json, inserts)
        section = {"Figure": figure, "Overview": overview_dict, "Content": html}

        return section
