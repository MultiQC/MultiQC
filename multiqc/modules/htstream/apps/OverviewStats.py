from multiqc.plots import linegraph


class OverviewStats:
    def read_and_basepair_reduction(self, json, app_list):
        line_config = {
            "id": "htstream_overview_reduction",
            "smooth_points_sumcounts": False,
            "categories": True,
            "title": "HTStream: Read and Basepair Reduction",
            "tt_label": "{point.y:.2}",
            "ylab": "Counts",
            "data_labels": [
                {"name": "Reads", "ylab": "Counts", "xlab": "Tool"},
                {"name": "Basepairs", "ylab": "Counts", "xlab": "Tool"},
            ],
        }

        data = [{}, {}]

        # Initialize lists for sample and app order
        samples = list(json["Pipeline Input"].keys())
        app_list = ["Pipeline Input"] + app_list  # preserves order of elements

        for samp in samples:
            # initilize read and bp line dicts
            data[0][samp] = {}
            data[1][samp] = {}

            # Iterate through app list, if desired app is found,
            #   grab total read counts and bp counts
            #   and add them to line graphs.

            for app in app_list:
                if app == "Pipeline Input":
                    io = "Input"

                else:
                    io = "Output"

                total_reads = json[app][samp][io + "_Reads"]
                data[0][samp][app] = total_reads

                total_bps = json[app][samp][io + "_Bps"]
                data[1][samp][app] = total_bps

        return linegraph.plot(data, line_config)

    def execute(self, json, app_list):
        return self.read_and_basepair_reduction(json, app_list)
