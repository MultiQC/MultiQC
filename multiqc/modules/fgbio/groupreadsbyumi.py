""" MultiQC module to parse output from fgbio GroupReadsByUmi
"""


import logging

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class GroupReadsByUmiMixin:
    def parse_groupreadsbyumi(self):
        # Parse data
        self.parse_groupreadsbyumi_log()

        # Make a plot if we have anything
        if len(self.fgbio_umi_data) > 0:
            self.parse_groupreadsbyumi_plot()

        # Return the number of samples parsed
        return len(self.fgbio_umi_data)

    def parse_groupreadsbyumi_log(self):
        umi_data = dict()
        umi_data_normed = dict()

        for f in self.find_log_files("fgbio/groupreadsbyumi"):
            # add file to data sources
            self.add_data_source(f)
            family_size = []
            for line in f["f"].splitlines():
                if not line.startswith("family_size"):
                    family_size.append(tuple(line.split("\t")))

            umi_data[f["s_name"]] = {int(s): int(d[1]) for s, d in enumerate(family_size, 1)}
            umi_data_normed[f["s_name"]] = {int(s): float(d[2]) * 100.0 for s, d in enumerate(family_size, 1)}

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        # Filter samples
        self.fgbio_umi_data = self.ignore_samples(umi_data)
        self.fgbio_umi_data_normed = self.ignore_samples(umi_data_normed)

        # Write data to file
        self.write_data_file(self.fgbio_umi_data, "fgbio_umi")
        self.write_data_file(self.fgbio_umi_data_normed, "fgbio_umi_normed")

    def parse_groupreadsbyumi_plot(self):
        config = {
            "id": "fgbio-groupreadsbyumi-plot",
            "title": "fgbio: Family size count",
            "ylab": "Number of UMIs",
            "xlab": "Reads supporting UMI",
            "xmax": 15,
            "data_labels": [
                {
                    "name": "Counts",
                    "ylab": "Number of UMIs",
                },
                {
                    "name": "Percentages",
                    "ylab": "Percentage of sample",
                    "ysuffix": "%",
                },
            ],
        }

        self.add_section(
            name="GroupReadsByUmi statistics",
            anchor="fgbio-groupreadsbyumi",
            description="""During `GroupReadsByUmi` processing, family size count data is generated,
                             showing number of UMIs represented by a certain number of reads.""",
            helptext="""
            This tool groups reads together that appear to have come from the same original molecule.
            Reads are grouped by template, and then templates are sorted by the 5' mapping positions
            of the reads from the template, used from earliest mapping position to latest.
            Reads that have the same end positions are then sub-grouped by UMI sequence.

            The histogram shows tag family size counts or percentages.
            """,
            plot=linegraph.plot([self.fgbio_umi_data, self.fgbio_umi_data_normed], config),
        )
