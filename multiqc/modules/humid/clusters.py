import logging
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_reports(self):
    # To store the summary data
    self.clusters = dict()

    # Parse the output files
    parse_log_files(self)

    # Remove filtered samples
    self.clusters = self.ignore_samples(self.clusters)

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    if self.clusters:
        self.write_data_file(self.clusters, "multiqc_humid_clusters")
        add_to_humid_section(self)


def parse_log_files(self):
    for f in self.find_log_files("humid/clusters"):
        # There is no sample name in the log, so we use the root of the
        # file as sample name (since the filename is always stats.dat
        s_name = self.clean_s_name(f["root"], f)

        # process the file content
        d = {}
        for line in f["f"].splitlines():
            cluster_size, count = line.strip("\n").split(" ")
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
        "hide_empty": False,
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
