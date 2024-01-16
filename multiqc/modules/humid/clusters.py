import logging
from multiqc.plots import linegraph

log = logging.getLogger(__name__)

def parse_reports(self):
    # To store the summary data
    self.clusters = dict()

    # Parse the output files
    parse_log_files(self)

    # Remove filtered samples
    self.cluster = self.ignore_samples(self.clusters)

    log.info(f"Found {len(self.clusters)} reports")
    self.write_data_file(self.clusters, "multiqc_humid_clusters")

    add_to_humid_section(self)


def parse_log_files(self):
    for f in self.find_log_files("humid/clusters"):
        # There is no sample name in the log, so we use the root of the
        # file as sample name (since the filename is always stats.dat
        s_name = self.clean_s_name(f["root"], f)

        # process the file content
        d = {}
        for line in f["contents_lines"]:
            cluster_size, count = line.strip('\n').split(' ')
            d[int(cluster_size)] = int(count)

        # Got this far, data must be good
        if s_name in self.clusters:
            log.debug("Duplicate sample name found! Overwriting: {s_name}")
        self.clusters[s_name] = d
        self.add_data_source(f, s_name)

def add_to_humid_section(self):
    # Figure configuration
    plot_config = {
        "id": "humid-clusters",
        "title": "HUMID: Cluster statistics",
        "ylab": "Number of clusters",
        "xlab": "Cluster size",
        "logswitch": True,
        "logswitch_active": True,
        "hide_zero_cats": False,
    }
    self.add_section(
        name="Cluster statistics",
        anchor="humid-cluster-section",
        description="""
            Cluster statistics per sample. For every cluster in the trie,
            the number of reads that fall into this cluster are shown.
            """,
        plot=linegraph.plot(self.counts, plot_config),
    )
