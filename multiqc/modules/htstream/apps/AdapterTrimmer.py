from collections import OrderedDict
import logging

from multiqc.plots import bargraph

#################################################

""" AdapterTrimmer submodule for HTStream charts and graphs """

#################################################


class AdapterTrimmer:
    ########################
    # Info about App
    def __init__(self):
        self.info = (
            "Trims adapters which are sequenced when the fragment insert length is shorter than the read length."
        )
        self.type = "bp_reducer"

    ########################
    # Table Function
    def bargraph(self, json, bps, index):
        # config dict for bar graph
        config = {
            "title": "HTStream: AdapterTrimmer Trimmed Basepairs Bargraph",
            "id": "htstream_adaptertrimmer_bargraph_" + index,
            "ylab": "Percentage of Total Basepairs",
            "cpswitch": False,
            "data_labels": [
                {"name": "Basepairs Lost", "ylab": "Percentage of Total Basepairs"},
                {"name": "Reads with Adapters", "ylab": "Percentage of Total Reads"},
                {"name": "Average Number of Basepairs Trimmed", "ylab": "Basepairs"},
            ],
        }

        # Header
        html = ""

        # returns nothing if no reads were trimmed.
        if bps == 0:
            html = '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from any sample. </div>'
            figure = None

        else:
            perc_data = {}
            read_data = {}
            bp_data = {}

            # Create dictionaries for multidataset bargraphs
            for key in json:
                perc_data[key] = {"Perc_bp_lost": json[key]["At_%_BP_Lost" + index]}
                read_data[key] = {"Perc_adapters": json[key]["At_%_Adapters" + index]}
                bp_data[key] = {"Avg_adapter": json[key]["At_Avg_BP_Trimmed" + index]}

            # Create categories for multidataset bargraph
            cats = [OrderedDict(), OrderedDict(), OrderedDict()]
            cats[0]["Perc_bp_lost"] = {"name": "Percentage"}
            cats[1]["Perc_adapters"] = {"name": "Percentage"}
            cats[2]["Avg_adapter"] = {"name": "Basepairs"}

            figure = bargraph.plot([perc_data, read_data, bp_data], cats, config)

        return figure, html

    ########################
    # Main Function
    def execute(self, json, index):
        stats_json = OrderedDict()
        overview_dict = {}
        total = 0

        for key in json.keys():
            frag_in = (2 * json[key]["Paired_end"]["in"]) + json[key]["Single_end"]["in"]
            bp_in = json[key]["Fragment"]["basepairs_in"]

            # try block to avoid zero division
            try:
                fract_bp_lost = (bp_in - json[key]["Fragment"]["basepairs_out"]) / bp_in
                perc_bp_lost = fract_bp_lost * 100

                fract_pe_bp_trimmed = (
                    json[key]["Paired_end"]["Read1"]["adapterBpTrim"]
                    + json[key]["Paired_end"]["Read2"]["adapterBpTrim"]
                ) / bp_in
                fract_pe_reads_trimmed = (
                    json[key]["Paired_end"]["Read1"]["adapterTrim"] + json[key]["Paired_end"]["Read2"]["adapterTrim"]
                ) / frag_in

                fract_se_bp_trimmed = json[key]["Single_end"]["adapterBpTrim"] / bp_in
                fract_se_reads_trimmed = json[key]["Single_end"]["adapterTrim"] / frag_in

            except ZeroDivisionError:
                fract_bp_lost = 0
                perc_bp_lost = 0

                fract_pe_bp_trimmed = 0
                fract_pe_reads_trimmed = 0

                fract_se_bp_trimmed = 0
                fract_se_reads_trimmed = 0

                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)

            # calculations for reads with adapters and bps trimmed
            adapter_reads = (
                json[key]["Single_end"]["adapterTrim"]
                + json[key]["Paired_end"]["Read1"]["adapterTrim"]
                + json[key]["Paired_end"]["Read2"]["adapterTrim"]
            )  # total reads trimmed
            bp_trimmed = (
                json[key]["Single_end"]["adapterBpTrim"]
                + json[key]["Paired_end"]["Read1"]["adapterBpTrim"]
                + json[key]["Paired_end"]["Read2"]["adapterBpTrim"]
            )  # total basepairs trimmed

            # if adapter trim is zero, so is the percentage and the avg basepair trimmed.
            #   This prevents division by zero error.
            if adapter_reads == 0 or frag_in == 0:
                perc_adapters = 0
                avg_bp_trimmed = 0

            else:
                perc_adapters = (adapter_reads / frag_in) * 100
                avg_bp_trimmed = bp_trimmed / frag_in

            # Accumulate avg bp trimmed
            total += avg_bp_trimmed

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Bp_Lost": fract_bp_lost,
                "Fraction_PE_Bp_Trimmed": fract_pe_bp_trimmed,
                "Fraction_PE_Read_Trimmed": fract_pe_reads_trimmed,
                "Fraction_SE_Bp_Trimmed": fract_se_bp_trimmed,
                "Fraction_SE_Read_Trimmed": fract_se_reads_trimmed,
            }

            # sample dictionary entry
            stats_json[key] = {
                "At_%_BP_Lost" + index: perc_bp_lost,
                "At_%_Adapters" + index: perc_adapters,
                "At_Avg_BP_Trimmed" + index: avg_bp_trimmed,
            }

        # sections and figure function calls
        figure, html = self.bargraph(stats_json, total, index)
        section = {"Figure": figure, "Overview": overview_dict, "Content": html}

        return section
