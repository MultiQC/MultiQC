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

        samples = dict()
        for f in list(log_files):
            if not f or "f" not in f:
                log.warning("Malformed log file object: {}".format(f))
            else:
                print("Parsing log file:", f["fn"])
                samples = self._parse_samples(f, samples)
                self.add_data_source(f)

        self._samples = samples

        # Filter to strip out ignored sample names
        self._samples = self.ignore_samples(self._samples)

        if len(self._samples) == 0:
            raise UserWarning

        self._add_table(samples)
        self._add_coverage_plot(samples)

        self.write_data_file(self._samples, "multiqc_mosaicatcher")
        log.info("Found {} samples".format(len(samples)))

    def _add_table(self, samples):
        self.add_section(
            plot=table.plot(
                data=samples,
                headers=self._setup_headers(),
                pconfig={"id": "mosaicatcher_table", "namespace": "MosaiCatcher"},
            )
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
            name="Average coverage per cell",
            anchor="mosaicatcher-coverage-id",
            description="Average coverage per cell",
            plot=coverage_plot,
        )

    def _setup_headers(self):
        unwanted_keys = ["medbin", "nb_p", "nb_a", "nb_r", "bam"]

        headers = OrderedDict()
        for k in list(self._samples.values())[0]:
            if k in unwanted_keys:
                continue
            headers[k] = {"title": k, "hidden": True}

        headers["sample"] = {
            "title": "Sample",
            "description": "Sample (has multiple cells)",
            "hidden": True,
            "namespace": "MosaiCatcher",
        }
        headers["cell"] = {
            "title": "Cell",
            "description": "Name of the cell.",
            "hidden": False,
            "namespace": "MosaiCatcher",
        }
        headers["mapped"] = {
            "title": "Mapped reads",
            "max": 4e6,
            "description": "Total number of reads seen",
            "format": "{:,d}",
            "scale": "RdYlGn",
            "hidden": False,
            "namespace": "MosaiCatcher",
        }
        headers["suppl"] = {
            "title": "Supplementary reads",
            "description": "Supplementary, secondary or QC-failed reads (filtered out)",
            "format": "{:,d}",
            "hidden": True,
            "namespace": "MosaiCatcher",
        }
        headers["dupl"] = {
            "title": "Duplicate reads",
            "max": 3e6,
            "description": "Reads filtered out as PCR duplicates",
            "format": "{:,d}",
            "scale": "OrRd",
            "hidden": False,
            "namespace": "MosaiCatcher",
        }
        headers["mapq"] = {
            "title": "Mapping quality",
            "description": "Reads filtered out due to low mapping quality",
            "format": "{:,d}",
            "hidden": True,
            "namespace": "MosaiCatcher",
        }
        headers["read2"] = {
            "title": "2nd read of pair",
            "description": "Reads filtered out as 2nd read of pair",
            "format": "{:,d}",
            "hidden": True,
            "namespace": "MosaiCatcher",
        }
        headers["good"] = {
            "title": "Good quality reads",
            "max": 8e5,
            "scale": "RdYlGn",
            "description": "Reads used for counting.",
            "format": "{:,d}",
            "hidden": False,
            "namespace": "MosaiCatcher",
        }
        headers["pass1"] = {
            "title": "Enough coverage?",
            "description": "Enough coverage? If false, ignore all columns from now",
            "hidden": False,
            "namespace": "MosaiCatcher",
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
                        value = "pass" if value == "1" else "fail"

                    # Add this column to the cell dictionary
                    samples[cell_name][key] = value

        return samples
