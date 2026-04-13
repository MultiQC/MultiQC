import csv
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        """MultiQC module for processing som.py output"""
        super(MultiqcModule, self).__init__(
            name="som.py",
            anchor="sompy",
            href="https://github.com/Illumina/hap.py/blob/master/doc/sompy.md",
            info=("Benchmarks somatic variant calls against gold standard truth datasets."),
            # No publication / DOI // doi=
        )

        self.add_software_version(None)

        self.sompy_raw_sample_names = set()
        self.sompy_combined_data = dict()
        self.sompy_indel_data = dict()
        self.sompy_snv_data = dict()

        for f in self.find_log_files("sompy"):
            self.parse_file(f)
            self.add_data_source(f)

        if len(self.sompy_raw_sample_names) == 0:
            raise ModuleNoSamplesFound

        log.info("Found %s sompy reports", len(self.sompy_raw_sample_names))

        helptext = (
            "No plots are generated, as som.py is generally run on single"
            " control samples (HD757, etc.). Ideally, precision, recall and"
            " F1 Score should all be as close to 1 as possible."
        )

        self.add_section(
            name="Combined",
            anchor="sompy-combined-plot",
            helptext=helptext,
            plot=table.plot(
                self.sompy_combined_data,
                self.generate_table_headers("_combined"),
                pconfig={
                    "id": "sompy_combined_plot",
                    "title": "som.py: combined",
                },
            ),
        )
        self.add_section(
            name="Indel",
            anchor="sompy-indel-plot",
            helptext=helptext,
            plot=table.plot(
                self.sompy_indel_data,
                self.generate_table_headers("_indel"),
                pconfig={
                    "id": "sompy_indel_plot",
                    "title": "som.py: indel",
                },
            ),
        )
        self.add_section(
            name="SNV",
            anchor="sompy-snv-plot",
            helptext=helptext,
            plot=table.plot(
                self.sompy_snv_data,
                self.generate_table_headers("_snv"),
                pconfig={
                    "id": "sompy_snv_plot",
                    "title": "som.py: SNV",
                },
            ),
        )

        self.write_data_file(self.sompy_combined_data, "multiqc_sompy_combined_data")
        self.write_data_file(self.sompy_indel_data, "multiqc_sompy_indel_data")
        self.write_data_file(self.sompy_snv_data, "multiqc_sompy_snv_data")

    def generate_table_headers(self, suffix: str = "") -> dict:
        """
        Generates the dict of table header metadata
        """
        header = {
            "unk": {
                "title": "Not assessed calls",
                "description": "Number of non-assessed query calls",
                "min": 0,
                "max": 1,
                "hidden": True,
                "format": "{:.4f}",
            },
            "total.truth": {
                "title": "Truth: Total",
                "description": "Total number of truth variants",
                "format": None,
                "hidden": True,
            },
            "total.query": {
                "title": "Query: Total",
                "description": "Total number of query calls",
                "format": None,
                "hidden": True,
            },
            "tp": {
                "title": "True Positive Variants",
                "description": "Number of true-positive calls",
                "scale": "Greens",
                "format": None,
            },
            "fn": {
                "title": "False Negative Variants",
                "description": "Calls in truth without matching query call",
                "scale": "Reds",
                "format": None,
            },
            "fp": {
                "title": "False Positive Variants",
                "description": "Number of false-positive calls",
                "format": None,
                "scale": "Reds",
                "hidden": True,
            },
            "recall": {
                "title": "Recall",
                "description": ("Recall for truth variant representation = TRUTH.TP / (TRUTH.TP + TRUTH.FN)"),
                "min": 0,
                "max": 1,
                "cond_formatting_rules": {
                    "verygreen": [
                        {"gte": 0.99},
                    ],
                    "green": [{"lt": 0.99}, {"gt": 0.98}],
                    "amber": [{"lt": 0.98}, {"gt": 0.90}],
                    "red": [{"lt": 0.90}],
                },
                "cond_formatting_colours": [
                    {"red": "#D2222D"},
                    {"amber": "#FFBF00"},
                    {"green": "#238823"},
                    {"verygreen": "#007000"},
                ],
                "format": "{:.4f}",
            },
            "precision": {
                "title": "Precision",
                "description": ("Precision of query variants = QUERY.TP / (QUERY.TP + QUERY.FP)"),
                "min": 0,
                "max": 1,
                "format": "{:.4f}",
                "hidden": True,
            },
        }

        return {f"{k}{suffix}": v for k, v in header.items()}

    def parse_file(self, f):
        """
        Reads the combined, indel and SNV data from sompy output file
        """
        if self.is_ignore_sample(f["s_name"]):
            return

        if f["s_name"] in self.sompy_raw_sample_names:
            log.warning(
                "Duplicate sample name found in %s! Overwriting: %s",
                f["root"],
                f["s_name"],
            )

        self.sompy_raw_sample_names.add(f["s_name"])

        reader = csv.DictReader(f["f"].split("\n"))

        for row in reader:
            row_id = f"{f['s_name']}_{row['type']}"

            if row["type"] == "records":
                if row_id not in self.sompy_combined_data:
                    self.sompy_combined_data[row_id] = {"sample_id": f["s_name"]}

                for fn in reader.fieldnames:
                    self.sompy_combined_data[row_id][fn + "_combined"] = row[fn]

            if row["type"] == "indels":
                if row_id not in self.sompy_indel_data:
                    self.sompy_indel_data[row_id] = {"sample_id": f["s_name"]}

                for fn in reader.fieldnames:
                    self.sompy_indel_data[row_id][fn + "_indel"] = row[fn]

            if row["type"] == "SNVs":
                if row_id not in self.sompy_snv_data:
                    self.sompy_snv_data[row_id] = {"sample_id": f["s_name"]}

                for fn in reader.fieldnames:
                    self.sompy_snv_data[row_id][fn + "_snv"] = row[fn]
