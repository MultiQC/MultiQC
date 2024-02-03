import logging
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_reports(self):
    # To store the summary data
    self.counts = dict()

    # Parse the output files
    parse_log_files(self)

    # Remove filtered samples
    self.counts = self.ignore_samples(self.counts)

    if self.counts:
        self.write_data_file(self.counts, "multiqc_humid_counts")
        add_to_humid_section(self)


def parse_log_files(self):
    for f in self.find_log_files("humid/counts"):
        # There is no sample name in the log, so we use the root of the
        # file as sample name (since the filename is always stats.dat
        s_name = self.clean_s_name(f["root"], f)

        # process the file content
        d = {}
        for line in f["contents_lines"]:
            nr_reads, count = line.strip("\n").split(" ")
            d[int(nr_reads)] = int(count)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Got this far, data must be good
        if s_name in self.counts:
            log.debug("Duplicate sample name found! Overwriting: {s_name}")
        self.counts[s_name] = d
        self.add_data_source(f, s_name)


def add_to_humid_section(self):
    # Figure configuration
    plot_config = {
        "id": "humid-counts",
        "title": "HUMID: Count statistics",
        "ylab": "Number of nodes",
        "xlab": "Number of identical reads in a node",
        "logswitch": True,
        "logswitch_active": True,
        "hide_zero_cats": False,
    }
    self.add_section(
        name="Counts statistics",
        anchor="humid-count-section",
        description="""
            Count statistics per sample. This shows the number of times every
            read has been seen (reads must always match exactly to be
            counted here.)
            """,
        plot=linegraph.plot(self.counts, plot_config),
    )
