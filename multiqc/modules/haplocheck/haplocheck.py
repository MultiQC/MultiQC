import csv
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Haplocheck",
            anchor="haplocheck",
            href="https://github.com/genepi/haplocheck",
            info=(
                "Haplocheck detects contamination patterns in mtDNA AND WGS sequencing studies by analyzing "
                "the mitochondrial DNA."
            ),
            doi="10.1101/gr.256545.119",
        )
        log_files = self.find_log_files("haplocheck", filehandles=True)
        samples = {}
        for f in log_files:
            samples = self._parse_samples(f, samples)
            self.add_data_source(f)

        self._samples = samples

        self._add_table(samples)
        self._add_coverage_plot(samples)
        self._add_heteroplasmy_level_plot(samples)

        self.write_data_file(self._samples, "multiqc_haplocheck")
        log.info("Found {} samples".format(len(samples)))

    def _add_table(self, samples):
        self.add_section(
            description=(
                "For further information regarding the output values, refer to the "
                "[tool's documentation](https://mitoverse.readthedocs.io/haplocheck/haplocheck/#textual-report-file)."
            ),
            plot=table.plot(data=samples, headers=self._setup_headers()),
        )

    def _add_coverage_plot(self, samples):
        coverage_plot = bargraph.plot(
            data=samples,
            cats=["Sample Coverage"],
            pconfig={
                "id": "haplocheck-coverage",
                "title": "Haplocheck: coverages",
                "xlab": "Sample",
                "ylab": "Average Coverage",
                "tt_suffix": "x",
            },
        )
        self.add_section(
            name="Average coverage per sample",
            anchor="haplocheck-coverage-id",
            description="Average coverage per sample",
            plot=coverage_plot,
        )

    def _add_heteroplasmy_level_plot(self, samples):
        coverage_plot = bargraph.plot(
            data=samples,
            cats=["Major Heteroplasmy Level", "Minor Heteroplasmy Level"],
            pconfig={
                "id": "haplocheck-het-level",
                "title": "Haplocheck: Heteroplasmy Level",
                "xlab": "Sample",
                "ylab": "Heteroplasmy Level",
                "tt_decimals": 2,
            },
        )
        self.add_section(
            name="Heteroplasmy level per sample",
            anchor="haplocheck-het-level-id",
            description="Heteroplasmy level per sample",
            plot=coverage_plot,
        )

    def _setup_headers(self):
        headers = OrderedDict()
        for k in list(self._samples.values())[0]:
            headers[k] = {"title": k, "hidden": True}

        headers["Contamination Status"] = {
            "title": "Contamination Status",
            "description": "This column can either be YES, NO or ND (not detectable)",
            "hidden": False,
            "bgcols": {"YES": "#f8d7da", "ND": "#fff3cd", "NO": "#d1e7dd"},
        }
        headers["Contamination Level"] = {
            "title": "Contamination Level",
            "description": "The detected contamination level for the sample",
            "bgcols": {"ND": "#d1e7dd"},
            "hidden": False,
        }
        headers["Major Heteroplasmy Level"] = {
            "title": "Major Heteroplasmy Level",
            "max": 1,
            "min": 0,
            "format": "{:,.3f}",
            "hidden": True,
        }

        headers["Minor Heteroplasmy Level"] = {
            "title": "Minor Heteroplasmy Level",
            "max": 1,
            "min": 0,
            "format": "{:,.3f}",
            "hidden": True,
        }
        headers["Major Haplogroup"] = {
            "title": "Major Haplogroup",
            "hidden": False,
        }
        headers["Minor Haplogroup"] = {
            "title": "Minor Haplogroup",
            "description": "Only consider it when contamination is detected.",
            "hidden": False,
        }
        return headers

    @staticmethod
    def _parse_samples(f, samples):
        for row in csv.DictReader(f["f"], delimiter="\t"):
            sample_name = row["Sample"]
            row.pop("Sample")
            if sample_name in samples:
                log.warning("Duplicate sample name found! Overwriting: {}".format(sample_name))
            samples[sample_name] = row

        return samples
