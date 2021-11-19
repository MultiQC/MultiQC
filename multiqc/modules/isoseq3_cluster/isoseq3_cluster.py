from multiqc.modules.base_module import BaseMultiqcModule
import logging
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph
import re
import pandas as pd

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="isoseq3 cluster",
            anchor="isoseq3_cluster",
            href="https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md#step-4---clustering",
            info=" is a Pacbio tool that generates Hifi reads CCS with a single linkage clustering method.",
            doi="",
        )

        self.isoseq3_cluster_data = dict()
        self.parse_log()

        self.isoseq3_cluster_data = self.ignore_samples(self.isoseq3_cluster_data)

        # If we found no data
        if not self.isoseq3_cluster_data:
            raise UserWarning

        log.info("Found {} reports".format(len(self.isoseq3_cluster_data)))

        self.write_data_file(self.isoseq3_cluster_data, "multiqc_isoseq3_cluster_report")
        self.add_general_stats()
        self.add_sections()

    def parse_log(self):
        for f in self.find_log_files("isoseq3_cluster/csv", filehandles=True):
            data = pd.read_csv(f["f"])
            filename = f["s_name"]
            self.isoseq3_cluster_data[filename] = data.value_counts("cluster_id")
            self.add_data_source(f)

    def add_general_stats(self):
        gstats_data = {}
        for s_name, attrs in self.isoseq3_cluster_data.items():
            gstats_data[s_name] = {}
            gstats_data[s_name]["n_cluster"] = len(attrs.to_dict())
            gstats_data[s_name]["mean_cluster_size"] = attrs.mean()

        headers = OrderedDict()
        headers["n_cluster"] = {
            "title": "Number of cluster",
            "description": "Number of clusters created during the clustering. (1 cluster = 1 transcript)",
            "format": "{:,}",
            "scale": "spectral",
        }
        headers["mean_cluster_size"] = {
            "title": "Mean cluster size",
            "description": "Avearge number of CCS clusterd to form transcripts.",
            "format": "{:,0f}",
            "scale": "RdYlGn",
        }

        self.general_stats_addcols(gstats_data, headers)

    def add_sections(self):
        self.add_section_cluster_size_bargraph()

    def add_section_cluster_size_bargraph(self):
        plot_data = dict()
        for s_name, data in self.isoseq3_cluster_data.items():
            df = data.to_frame()
            df.columns = ["n_CCS"]
            df = df.value_counts("n_CCS").to_frame()
            df.columns = ["count"]
            df.reset_index(inplace=True)

            plot_data[s_name] = dict()
            plot_data[s_name]["=2"] = df[df["n_CCS"] == 2]["count"][0]
            plot_data[s_name]["3-10"] = sum(df[(df["n_CCS"] > 2) & (df["n_CCS"] < 11)]["count"])
            plot_data[s_name]["11-100"] = sum(df[(df["n_CCS"] > 10) & (df["n_CCS"] < 101)]["count"])
            plot_data[s_name][">100"] = sum(df[df["n_CCS"] > 100]["count"])

        cats = OrderedDict()
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
            "id": "isoseq3_cluster_cluster_size_histogram",  # HTML ID used for plot
            "title": "isoseq3 cluster: Histogram of cluster size.",  # Plot title - should be in format "Module Name: Plot Title"
            "xlab": "Cluster size",  # X axis label
            "ylab": "Count",  # Y axis label
        }

        self.add_section(
            name="Cluster size distribution",
            anchor="cluster-distribution",
            description="A distribution of cluster size (number of CC clustered to form one Hifi read).",
            helptext="""
                        The csv file produced by isoseq3 cluster shows which CCS have been clustered together to form \\  
                        one Hifi reads. The bargraph represent the distribution of the cluster size using four \\
                        categories : -2, 3-10, 11-100, >100.
                        """,
            plot=bargraph.plot(plot_data, cats, config),
        )
