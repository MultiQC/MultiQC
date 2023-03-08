import csv
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Haplocheck",
            anchor="haplocheck",
            href="https://github.com/genepi/haplocheck",
            info="Haplocheck detects contamination patterns in mtDNA AND WGS sequencing studies by analyzing the mitochondrial DNA.",
            doi="10.1101/gr.256545.119",
        )
        log_files = self.find_log_files("haplocheck", filehandles=True)
        samples = {}
        for f in log_files:
            samples = self._parse_samples(f, samples)
            self.add_data_source(f)

        self._samples = samples
        self.add_section(plot=table.plot(data=self._samples, headers=self._setup_headers()))
        self.write_data_file(self._samples, "multiqc_haplocheck")

    def _setup_headers(self):
        headers = OrderedDict()
        for k in list(self._samples.values())[0]:
            headers[k] = {"title": k, "hidden": True}

        headers["Contamination Status"] = {
            "title": "Contamination Status",
            "hidden": False,
            "bgcols": {"YES": "#f8d7da", "ND": "#fff3cd", "NO": "#d1e7dd"},
        }
        headers["Contamination Level"] = {
            "title": "Contamination Level",
            "hidden": False,
        }
        headers["Major Heteroplasmy Level"] = {
            "title": "Major Heteroplasmy Level",
            "max": 1,
            "min": 0,
            "format": "{:,.3f}",
            "hidden": False,
        }

        headers["Minor Heteroplasmy Level"] = {
            "title": "Minor Heteroplasmy Level",
            "max": 1,
            "min": 0,
            "format": "{:,.3f}",
            "hidden": False,
        }
        headers["Major Haplogroup"] = {
            "title": "Major Haplogroup",
            "hidden": False,
        }
        headers["Minor Haplogroup"] = {
            "title": "Minor Haplogroup",
            "hidden": False,
        }
        return headers

    @staticmethod
    def _parse_samples(f, samples):
        for row in csv.DictReader(f["f"], delimiter="\t"):
            sample_name = row["Sample"].replace(".raw", "")
            row.pop("Sample")
            if sample_name in samples:
                log.debug("Duplicate sample name found! Overwriting: {}".format(sample_name))
            samples[sample_name] = row

        return samples
