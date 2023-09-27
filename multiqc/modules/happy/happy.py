"""MultiQC module to parse output from OUS variant calling pipeline"""


import csv
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

log = logging.getLogger(__name__)
VAR_TYPES = ("INDEL", "SNP")
VAR_FILTERS = ("ALL", "PASS")
METRICS = OrderedDict(
    [
        ("Recall", "Recall for truth variant representation = TRUTH.TP / (TRUTH.TP + TRUTH.FN)"),
        ("Precision", "Precision of query variants = QUERY.TP / (QUERY.TP + QUERY.FP)"),
        (
            "F1_Score",
            "Harmonic mean of precision and recall = 2METRIC.RecallMetric.Precision/(METRIC.Recall + METRIC.Precision)",
        ),
        ("Frac_NA", "Fraction of non-assessed query calls = QUERY.UNK / QUERY.TOTAL"),
    ]
)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        """MultiQC module for processing hap.py output logs"""
        super(MultiqcModule, self).__init__(
            name="hap.py",
            anchor="happy",
            href="https://github.com/Illumina/hap.py",
            info="is a set of programs based on htslib to benchmark variant calls against gold standard truth datasets.",
            # No publication / DOI // doi=
        )

        self.happy_raw_sample_names = set()
        self.happy_indel_data = dict()
        self.happy_snp_data = dict()

        for f in self.find_log_files("happy", filehandles=True):
            self.parse_file(f)
            self.add_data_source(f)

        if len(self.happy_raw_sample_names) == 0:
            raise UserWarning

        # print number of happy reports found and parsed
        log.info("Found {} reports".format(len(self.happy_raw_sample_names)))

        # Write parsed report data to file
        self.write_data_file(self.happy_indel_data, "multiqc_happy_indel_data", data_format="json")
        self.write_data_file(self.happy_snp_data, "multiqc_happy_snp_data", data_format="json")

        # add sections with the values from the indel and snp happy output
        self.add_section(
            name="INDEL",
            anchor="happy-indel-plot",
            description="The default shown fields should give the best overview of quality, but there are many other hidden fields available.",
            helptext="""
                No plots are generated, as hap.py is generally run on single control samples (NA12878, etc.)

                Ideally, precision, recall and F1 Score should all be as close to 1 as possible.
            """,
            plot=table.plot(self.happy_indel_data, self.gen_headers("_indel")),
        )

        self.add_section(
            name="SNP",
            anchor="happy-snp-plot",
            description="The default shown fields should give the best overview of quality, but there are many other hidden fields available.",
            helptext="""
                No plots are generated, as hap.py is generally run on single control samples (NA12878, etc.)

                Ideally, precision, recall and F1 Score should all be as close to 1 as possible.
            """,
            plot=table.plot(self.happy_snp_data, self.gen_headers("_snp")),
        )

    def parse_file(self, f):
        # Check that we're not ignoring this sample name
        if self.is_ignore_sample(f["s_name"]):
            return

        if f["s_name"] in self.happy_raw_sample_names:
            log.warning("Duplicate sample name found in {}! Overwriting: {}".format(f["root"], f["s_name"]))
        self.happy_raw_sample_names.add(f["s_name"])

        rdr = csv.DictReader(f["f"])
        for row in rdr:
            row_id = "{}_{}_{}".format(f["s_name"], row["Type"], row["Filter"])
            if row["Type"] == "INDEL":
                if row_id not in self.happy_indel_data:
                    self.happy_indel_data[row_id] = {"sample_id": f["s_name"]}
                # add suffix for headers to differentiate between indel and snp
                for fn in rdr.fieldnames:
                    self.happy_indel_data[row_id][fn + "_indel"] = row[fn]

            elif row["Type"] == "SNP":
                if row_id not in self.happy_snp_data:
                    self.happy_snp_data[row_id] = {"sample_id": f["s_name"]}
                for fn in rdr.fieldnames:
                    self.happy_snp_data[row_id][fn + "_snp"] = row[fn]

    def gen_headers(self, suffix=""):
        h = OrderedDict()
        h["METRIC.Recall"] = {  # string must match headers in the input file
            "title": "Recall",  # whatever string to be displayed in the html report table
            "description": "Recall for truth variant representation = TRUTH.TP / (TRUTH.TP + TRUTH.FN)",
            "min": 0,
            "max": 1,
            "format": "{:.4f}",
        }

        h["METRIC.Precision"] = {
            "title": "Precision",
            "description": "Precision of query variants = QUERY.TP / (QUERY.TP + QUERY.FP)",
            "min": 0,
            "max": 1,
            "format": "{:.4f}",
        }

        h["METRIC.Frac_NA"] = {
            "title": "Fraction NA",
            "description": "Fraction of non-assessed query calls = QUERY.UNK / QUERY.TOTAL",
            "min": 0,
            "max": 1,
            "format": "{:.4f}",
        }

        h["METRIC.F1_Score"] = {
            "title": "F1 Score",
            "description": "Harmonic mean of precision and recall = 2METRIC.RecallMetric.Precision/(METRIC.Recall + METRIC.Precision)",
            "min": 0,
            "max": 1,
            "format": "{:.4f}",
        }

        h["TRUTH.TOTAL"] = {
            "title": "Truth: Total",
            "description": "Total number of truth variants",
            "format": None,
            "hidden": True,
        }

        h["TRUTH.TP"] = {
            "title": "Truth: True Positive",
            "description": "Number of true-positive calls in truth representation (counted via the truth sample column)",
            "format": None,
            "hidden": True,
        }

        h["TRUTH.FN"] = {
            "title": "Truth: False Negative",
            "description": "Number of false-negative calls = calls in truth without matching query call",
            "format": None,
            "hidden": True,
        }

        h["QUERY.TOTAL"] = {
            "title": "Query: Total",
            "description": "Total number of query calls",
            "format": None,
            "hidden": True,
        }

        h["QUERY.TP"] = {
            "title": "Query: True Positive",
            "description": "Number of true positive calls in query representation (counted via the query sample column)",
            "format": None,
            "hidden": True,
        }

        h["QUERY.FP"] = {
            "title": "Query: False Positive",
            "description": "Number of false-positive calls in the query file (mismatched query calls within the confident regions)",
            "format": None,
            "hidden": True,
        }

        h["QUERY.UNK"] = {
            "title": "Query: Unknown",
            "description": "Number of query calls outside the confident regions",
            "format": None,
            "hidden": True,
        }

        h["FP.gt"] = {
            "title": "False Positive genotype",
            "description": "Number of genotype mismatches (alleles match, but different zygosity)",
            "hidden": True,
        }

        h["FP.al"] = {
            "title": "False Positive allele",
            "description": "Number of allele mismatches (variants matched by position and not by haplotype)",
            "hidden": True,
        }

        h["TRUTH.TOTAL.TiTv_ratio"] = {
            "title": "Truth: Total TiTv ratio",
            "description": "Transition / Transversion ratio for all truth variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["TRUTH.TOTAL.het_hom_ratio"] = {
            "title": "Truth: Total het/hom ratio",
            "description": "Het/Hom ratio for all truth variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["TRUTH.FN.TiTv_ratio"] = {
            "title": "Truth: False Negative TiTv ratio",
            "description": "Transition / Transversion ratio for false-negative variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["TRUTH.FN.het_hom_ratio"] = {
            "title": "Truth: False Negative het/hom ratio",
            "description": "Het/Hom ratio for false-negative variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["TRUTH.TP.TiTv_ratio"] = {
            "title": "Truth: True Positive TiTv ratio",
            "description": "Transition / Transversion ratio for true positive variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["TRUTH.TP.het_hom_ratio"] = {
            "title": "Truth: True Positive het/hom ratio",
            "description": "Het/Hom ratio for true positive variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["QUERY.FP.TiTv_ratio"] = {
            "title": "Query: False Positive TiTv ratio",
            "description": "Transition / Transversion ratio for false positive variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["QUERY.FP.het_hom_ratio"] = {
            "title": "Query: False Positive het/hom ratio",
            "description": "Het/Hom ratio for false-positive variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["QUERY.TOTAL.TiTv_ratio"] = {
            "title": "Query: total TiTv ratio",
            "description": "Transition / Transversion ratio for all query variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["QUERY.TOTAL.het_hom_ratio"] = {
            "title": "Query: Total het/hom ratio",
            "description": "Het/Hom ratio for all query variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["QUERY.TP.TiTv_ratio"] = {
            "title": "Query: True Positive TiTv ratio",
            "description": "Transition / Transversion ratio for true positive variants (query representation)",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["QUERY.TP.het_hom_ratio"] = {
            "title": "Query: True Positive het/hom ratio",
            "description": "Het/Hom ratio for true positive variants (query representation)",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["QUERY.UNK.TiTv_ratio"] = {
            "title": "Query: Uknown TiTv ratio",
            "description": "Transition / Transversion ratio for unknown variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["QUERY.UNK.het_hom_ratio"] = {
            "title": "Query: Unknown het/hom ratio",
            "description": "Het/Hom ratio for unknown variants",
            "hidden": True,
            "format": "{:.4f}",
        }

        h["Subset.Size"] = {
            "title": "Subset Size",
            "description": "the number of nucleotides contained in the current subset",
            "format": None,
            "hidden": True,
        }

        h["Subset.IS_CONF.Size"] = {
            "title": "Subset Confident Size",
            "description": "This gives the number of confident bases (-f regions) in the current subset",
            "format": None,
            "hidden": True,
        }
        # rename column headers with '_indel' or '_snp' suffix
        headers = [k + suffix for k in h.keys()]
        # recreate the ordered dictionary with all headers and information
        header_dict = OrderedDict(zip(headers, h.values()))
        return header_dict
