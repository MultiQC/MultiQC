import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="MosaiCatcher",
            anchor="mosaicatcher",
            href="https://github.com/friendsofstrandseq/mosaicatcher",
            info=(
                "Counts strand-seq reads and classifies strand states "
                "of each chromosome in each cell using a Hidden Markov Model."
            ),
            doi="10.1038/s41587-019-0366-x",
        )

        data_by_sample = dict()
        for f in list(self.find_log_files("mosaicatcher", filehandles=False)):
            data_by_sample = self._parse_samples(f, data_by_sample)
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self._add_table(data_by_sample)
        self._add_coverage_plot(data_by_sample)

        self.write_data_file(data_by_sample, "multiqc_mosaicatcher")

    def _add_table(self, data_by_sample):
        headers = {}
        headers["sample"] = {
            "title": "Sample",
            "description": "Sample (has multiple cells)",
            "hidden": True,
        }
        headers["cell"] = {
            "title": "Cell",
            "description": "Name of the cell.",
            "hidden": False,
        }
        headers["mapped"] = {
            "title": "Mapped",
            "description": "Total number of reads seen",
            "format": "{:,d}",
            "scale": "RdYlGn",
            "shared_key": "reads",
        }
        headers["suppl"] = {
            "title": "Supplementary",
            "description": "Supplementary, secondary or QC-failed reads (filtered out)",
            "format": "{:,d}",
            "hidden": True,
            "shared_key": "reads",
        }
        headers["dupl"] = {
            "title": "Duplicate",
            "description": "Reads filtered out as PCR duplicates",
            "format": "{:,d}",
            "scale": "OrRd",
            "hidden": False,
            "shared_key": "reads",
        }
        headers["mapq"] = {
            "title": "Low MQ",
            "description": "Reads filtered out due to low mapping quality",
            "format": "{:,d}",
            "hidden": True,
            "shared_key": "reads",
        }
        headers["read2"] = {
            "title": "2nd in pair",
            "description": "Reads filtered out as 2nd read of pair",
            "format": "{:,d}",
            "hidden": True,
            "shared_key": "reads",
        }
        headers["good"] = {
            "title": "Good quality reads",
            "scale": "RdYlGn",
            "description": "Reads used for counting.",
            "format": "{:,d}",
            "shared_key": "reads",
        }
        headers["pass1"] = {
            "title": "Enough coverage?",
            "description": "Enough coverage? If false, ignore all columns from now",
            "hidden": False,
        }
        self.add_section(
            name="MosaiCatcher statistics",
            anchor="mosaicatcher-table-section",
            plot=table.plot(
                data=data_by_sample,
                headers=headers,
                pconfig={"id": "mosaicatcher-table", "title": "MosaiCatcher: statistics"},
            ),
        )

    def _add_coverage_plot(self, data_by_sample):
        self.add_section(
            name="Average coverage per cell",
            anchor="mosaicatcher-coverage-section",
            description="Average coverage per cell",
            plot=bargraph.plot(
                data=data_by_sample,
                cats=[
                    {
                        "mapq": {"name": "Low MQ"},
                        "dupl": {"name": "Duplicate"},
                        "good": {"name": "Good quality reads"},
                        "reads2": {"name": "2nd in pair"},
                        "suppl": {"name": "Supplementary"},
                    }
                ],
                pconfig={
                    "id": "mosaicatcher-coverage",
                    "title": "MosaiCatcher: coverage",
                    "height": 1024,
                    "ymax": 4e6,
                    "ylab": "Cells",
                },
            ),
        )

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
