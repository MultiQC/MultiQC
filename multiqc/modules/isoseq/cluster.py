"""
MultiQC submodule to parse output from Iso-Seq refine logs
"""

import csv
import logging
from collections import defaultdict

from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class ClusterMixin:
    def parse_cluster_log(self):
        data = dict()

        for f in self.find_log_files("isoseq/cluster-csv", filehandles=True):
            d = defaultdict(int)
            reader = csv.DictReader(f["f"])
            for row in reader:
                cluster_id = row.get("cluster_id", None)
                if cluster_id is not None:
                    d[cluster_id] += 1
            if d:
                data[f["s_name"]] = d
                self.add_data_source(f)
        self.write_data_file(data, "multiqc_isoseq_cluster_report")
        return data

    def add_general_stats_cluster(self, data):
        gstats_data = {}
        for s_name, attrs in data.items():
            gstats_data[s_name] = {}
            gstats_data[s_name]["n_cluster"] = len(attrs)
            gstats_data[s_name]["mean_cluster_size"] = sum(attrs.values()) / len(attrs)

        headers = dict()
        headers["n_cluster"] = {
            "title": "# clusters",
            "description": "Number of clusters created during the clustering. (1 cluster = 1 transcript)",
            "scale": "spectral",
        }
        headers["mean_cluster_size"] = {
            "title": "Mean cluster size",
            "description": "Average number of CCS clustered to form transcripts.",
            "format": "{:,0f}",
            "scale": "RdYlGn",
        }

        self.general_stats_addcols(gstats_data, headers)

    def add_section_cluster_size_bargraph(self, data):
        plot_data = dict()

        for s_name, data_dict in data.items():
            plot_data[s_name] = {"=2": 0, "3-10": 0, "11-100": 0, ">100": 0}
            value_counts = defaultdict(int)

            # Calculating value counts similar to df.value_counts("n_CCS") in pandas
            for key, count in data_dict.items():
                value_counts[count] += 1

            for n_CCS, count in value_counts.items():
                if n_CCS == 2:
                    plot_data[s_name]["=2"] = count
                elif 2 < n_CCS < 11:
                    plot_data[s_name]["3-10"] += count
                elif 10 < n_CCS < 101:
                    plot_data[s_name]["11-100"] += count
                elif n_CCS > 100:
                    plot_data[s_name][">100"] += count

        cats = dict()
        cats["=2"] = {
            "name": "2 CCS",
        }
        cats["3-10"] = {
            "name": "3 to 10 CCS",
        }
        cats["11-100"] = {
            "name": "11 to 100 CCS",
        }
        cats[">100"] = {
            "name": "More than 100 CCS",
        }

        config = {
            "id": "isoseq_cluster_cluster_size_histogram",  # HTML ID used for plot
            "title": "Iso-Seq cluster: Histogram of cluster size.",  # Plot title - should be in format "Module Name: Plot Title"
            "xlab": "Cluster size",  # X axis label
            "ylab": "Count",  # Y axis label
        }

        self.add_section(
            name="Cluster size distribution",
            anchor="cluster-distribution",
            description="A distribution of cluster size (number of CC clustered to form one Hifi read).",
            helptext="""
            The csv file produced by Iso-Seq cluster shows which CCS have been clustered together to form 
            one Hifi reads. The bargraph represent the distribution of the cluster size using four 
            categories : -2, 3-10, 11-100, >100.
            """,
            plot=bargraph.plot(plot_data, cats, config),
        )
