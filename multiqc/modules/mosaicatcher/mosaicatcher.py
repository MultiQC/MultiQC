import csv
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
        log_files = self.find_log_files("mosaicatcher", filehandles=True)
        samples = {}
        for f in log_files:
            samples = self._parse_samples(f, samples)
            self.add_data_source(f)

        self._samples = samples

        # self._add_table(samples)
        # self._add_coverage_plot(samples)
        # self._add_heteroplasmy_level_plot(samples)

        self.write_data_file(self._samples, "multiqc_mosaicatcher")
        log.info("Found {} samples".format(len(samples)))

    @staticmethod
    def _parse_samples(f, samples):
        for row in csv.DictReader(f["f"], delimiter="\t"):
            sample_name = row["Sample"]
            row.pop("Sample")
            if sample_name in samples:
                log.warning("Duplicate sample name found! Overwriting: {}".format(sample_name))
            samples[sample_name] = row

        return samples
