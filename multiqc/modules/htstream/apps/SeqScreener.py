from collections import OrderedDict
import logging

from multiqc.plots import bargraph

#################################################

""" SeqScreener submodule for HTStream charts and graphs """

#################################################


class SeqScreener:
    ########################
    # Info about App
    def __init__(self):
        self.info = "A simple sequence screening tool which uses a kmer lookup approach to identify reads from an unwanted source."
        self.type = "read_reducer"

    # Bargraph Function
    def bargraph(self, json, reads_screened, index):
        # config dict for bar graph
        config = {
            "title": "HTStream: Identified Reads Bargraph",
            "id": "htstream_seqscreener_bargraph_" + index,
            "cpswitch": False,
            "ylab": "Percentage of Total Reads",
            "data_labels": [
                {"name": "Percentage of Total", "ylab": "Percentage of Total Reads"},
                {"name": "Raw Counts", "ylab": "Reads"},
            ],
        }

        html = ""

        # returns nothing if no reads were trimmed.
        if reads_screened == 0:
            html += '<div class="alert alert-info"> <strong>Notice:</strong> No sequences were identified in any sample. </div>'
            figure = None

        else:
            perc_data = {}
            read_data = {}

            # Construct data for multidataset bargraph
            for key in json:
                perc_data[key] = {"Perc_PE": json[key]["Ss_PE_%_hits"], "Perc_SE": json[key]["Ss_SE_%_hits"]}
                read_data[key] = {"Reads_PE": json[key]["Ss_PE_hits"], "Reads_SE": json[key]["Ss_SE_hits"]}

            # Create categories for multidataset bargraph
            cats = [OrderedDict(), OrderedDict()]
            cats[0]["Perc_PE"] = {"name": "Paired End"}
            cats[0]["Perc_SE"] = {"name": "Single End"}
            cats[1]["Reads_PE"] = {"name": "Paired End"}
            cats[1]["Reads_SE"] = {"name": "Single End"}

            figure = bargraph.plot([perc_data, read_data], cats, config)

        return figure, html

    ########################
    # Main Function
    def execute(self, json, index):
        stats_json = OrderedDict()
        overview_dict = {}

        reads_screened = 0

        for key in json.keys():
            pe_hits = json[key]["Paired_end"]["hits"]
            se_hits = json[key]["Single_end"]["hits"]

            # Will fail if no PE data
            try:
                perc_pe_hits = (pe_hits / json[key]["Paired_end"]["in"]) * 100

            except ZeroDivisionError:
                perc_pe_hits = 0

            # Will fail if no SE data
            try:
                perc_se_hits = (se_hits / json[key]["Single_end"]["in"]) * 100

            except ZeroDivisionError:
                perc_se_hits = 0

            # for overview section
            try:
                fract_reads_lost = (json[key]["Fragment"]["in"] - json[key]["Fragment"]["out"]) / json[key]["Fragment"][
                    "in"
                ]
                fract_hits = (pe_hits + se_hits) / json[key]["Fragment"]["in"]

            except ZeroDivisionError:
                fract_reads_lost = 0
                fract_hits = 0

                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)

            reads_screened += pe_hits + se_hits

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Reads_Lost": fract_reads_lost,
                "Fraction_Hits": fract_hits,
            }

            # sample entry for stats dictionary
            stats_json[key] = {
                "Ss_PE_%_hits": perc_pe_hits,
                "Ss_PE_hits": pe_hits,
                "Ss_SE_%_hits": perc_se_hits,
                "Ss_SE_hits": se_hits,
            }

        # sections and figure function calls
        figure, html = self.bargraph(stats_json, reads_screened, index)
        section = {"Figure": figure, "Overview": overview_dict, "Content": html}

        return section
