import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="MosaiCatcher",
            anchor="mosaicatcher",
            href="https://github.com/friendsofstrandseq/mosaicatcher",
            info=(
                "Mosaicatcher counts Strand-seq reads and classifies strand states "
                "of each chromosome in each cell using a Hidden Markov Model."
            ),
            doi="10.1038/s41587-019-0366-x",
        )
        log_files = self.find_log_files("mosaicatcher", filehandles=False)
        if not log_files:
            log.warning("No log files found for MosaiCatcher module")
        else:
            print("Log files found:", log_files)

        samples = dict()
        for f in list(log_files):
            if not f or "f" not in f:
                log.warning("Malformed log file object: {}".format(f))
            else:
                print("Parsing log file:", f["fn"])
                samples = self._parse_samples(f, samples)
                self.add_data_source(f)

        self._samples = samples

        self._add_table(samples)
        self._add_coverage_plot(samples)

        self.write_data_file(self._samples, "multiqc_mosaicatcher")
        log.info("Found {} samples".format(len(samples)))

    def _add_table(self, samples):
        self.add_section(
            plot=table.plot(data=samples, headers=self._setup_headers()),
        )

    def _add_coverage_plot(self, samples):
        coverage_plot = bargraph.plot(
            data=samples,
            cats=["mapq", "dupl", "good", "reads2", "suppl"],
            pconfig={
                "id": "mosaicatcher-coverage",
                "title": "MosaiCatcher: coverage",
                "height": 1024,
                "ymax": 4e6,
                "ylab": "Cells",
            },
        )
        self.add_section(
            name="Average coverage per sample",
            anchor="mosaicatcher-coverage-id",
            description="Average coverage per sample",
            plot=coverage_plot,
        )

    def _setup_headers(self):
        headers = OrderedDict()
        for k in list(self._samples.values())[0]:
            headers[k] = {"title": k, "hidden": True}

        headers["mosaicatcher-sample"] = {
            "title": "sample",
            "description": "Sample (has multiple cells)",
            "hidden": True,
        }
        headers["mosaicatcher-cell"] = {
            "title": "cell",
            "description": "Name of the cell.",
            "hidden": False,
        }
        headers["mosaicatcher-mapped"] = {
            "title": "mapped",
            "max": 4e6,
            "description": "Total number of reads seen",
            "scale": "RdYlGn",
            "hidden": False,
        }
        headers["mosaicatcher-suppl"] = {
            "title": "suppl",
            "description": "Supplementary, secondary or QC-failed reads (filtered out)",
            "hidden": True,
        }
        headers["mosaicatcher-dupl"] = {
            "title": "dupl",
            "max": 3e6,
            "description": "Reads filtered out as PCR duplicates",
            "scale": "OrRd",
            "hidden": False,
        }
        headers["mosaicatcher-mapq"] = {
            "title": "mapq",
            "description": "Reads filtered out due to low mapping quality",
            "hidden": True,
        }
        headers["mosaicatcher-read2"] = {
            "title": "read2",
            "description": "Reads filtered out as 2nd read of pair",
            "hidden": True,
        }
        headers["mosaicatcher-good"] = {
            "title": "good",
            "max": 8e5,
            "scale": "RdYlGn",
            "description": "Reads used for counting.",
            "hidden": False,
        }
        headers["mosaicatcher-pass1"] = {
            "title": "pass1",
            "description": "Enough coverage? If false, ignore all columns from now",
            "hidden": False,
        }
        headers["mosaicatcher-nb_p"] = {
            "title": "nb_p",
            "description": "Negative Binomial parameter p. Constant for one sample.",
            "hidden": True,
        }
        headers["mosaicatcher-nb_r"] = {
            "title": "nb_r",
            "description": "Negative Binomial parameter r. We use NB(p,r/2) * NB(p,r/2) in WC states, but NB(p,(1-a)*r)*NB(p,a*r) in WW or CC states.",
            "hidden": True,
        }
        headers["mosaicatcher-nb_a"] = {
            "title": "nb_a",
            "description": "Negative Binomial parameter a (alpha) used for zero expectation (see above).",
            "hidden": True,
        }
        headers["mosaicatcher-bam"] = {
            "title": "bam",
            "description": "Bam file of this cell",
            "hidden": True,
        }
        return headers

    @staticmethod
    def _parse_samples(f, samples):
        # Create a dictionary to store the data

        for row in f["f"].split("\n"):
            tmp_row = row.split("\t")
            if len(tmp_row) != 14:
                continue

            if tmp_row[0] == "sample":
                header = tmp_row

            else:
                # Extract the cell name from the row
                cell_name = tmp_row[1]

                # Create a dictionary for this cell if it doesn't exist yet
                if cell_name not in samples:
                    samples[cell_name] = {}

                # Loop over the columns in the row and add them to the cell dictionary
                for i, value in enumerate(tmp_row):
                    # Use the header value as the key for this column
                    key = header[i]

                    # Skip the cell column since it's being used as the top-level key
                    if key == "cell":
                        continue

                    if key == "pass1":
                        value = True if value == "1" else False

                    # Add this column to the cell dictionary
                    samples[cell_name][key] = value

        return samples
