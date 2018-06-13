#!/usr/bin/env python

"""MultiQC module to parse output from OUS variant calling pipeline"""

from __future__ import print_function
from collections import OrderedDict
import json
import csv
import logging
import re

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)
VAR_TYPES = ("INDEL", "SNP")
VAR_FILTERS = ("ALL", "PASS")
METRICS = OrderedDict([
    ("Recall", "Recall for truth variant representation = TRUTH.TP / (TRUTH.TP + TRUTH.FN)"),
    ("Precision", "Precision of query variants = QUERY.TP / (QUERY.TP + QUERY.FP)"),
    ("F1_Score", "Harmonic mean of precision and recall = 2METRIC.RecallMetric.Precision/(METRIC.Recall + METRIC.Precision)"),
    ("Frac_NA", "Fraction of non-assessed query calls = QUERY.UNK / QUERY.TOTAL")
])


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):
        """ MultiQC module for processing hap.py output logs """
        super(MultiqcModule, self).__init__(
            name='hap.py',
            anchor='happy',
            href='https://github.com/Illumina/hap.py',
            info=' is a set of programs based on htslib to benchmark variant calls against gold standard truth datasets.'
        )

        self.general_stats_header = OrderedDict()
        self.general_stats_data = dict()
        self.happy_seen = set()
        self.happy_data = dict()

        for vtype in VAR_TYPES:
            for vfilt in VAR_FILTERS:
                for metric in METRICS:
                    hname = "{}/{}_METRIC.{}".format(vtype, vfilt, metric)
                    htitle = "{}/{} {}".format(vtype, vfilt, metric.replace('_', ' '))
                    hdesc = "{}/{} - {}".format(vtype, vfilt, METRICS[metric])
                    self.general_stats_header[hname] = {
                        "title": htitle,
                        "hidden": vfilt == "ALL",
                        "format": "{:.4f}",
                        "min": 0,
                        "max": 1,
                        "description": hdesc
                    }

        n_files = 0
        for f in self.find_log_files("happy", filehandles=True):
            f['s_name'] = self.clean_s_name(f['s_name'], f['root'])

            n = self.parse_log(f)
            if n > 0:
                n_files += 1
                self.add_data_source(f)
        log.info("Found {} hap.py log file(s)".format(n_files))
        if n_files == 0:
            raise UserWarning

        self.general_stats_addcols(self.general_stats_data, self.general_stats_header)

        if len(self.happy_data) > 0:
            self.write_data_file(self.happy_data, 'multiqc_happy_data', data_format="json")

            hp_headers = gen_headers()
            hp_plot = table.plot(self.happy_data, hp_headers)
            self.add_section(
                name="hap.py",
                anchor="happy-plot",
                plot=hp_plot
            )

    def parse_log(self, log_file):
        n = 0
        log.debug("Processing file {}".format(log_file["root"]))
        samp = log_file['s_name']

        if samp in self.happy_seen:
            log.warn("Duplicate sample name found in {}! Overwriting: {}".format(log_file['root'], log_file['s_name']))
        else:
            self.happy_seen.add(samp)
            self.general_stats_data[samp] = dict()

        rdr = csv.DictReader(log_file['f'])
        for row in rdr:
            row_id = "{}_{}_{}".format(samp, row["Type"], row["Filter"])
            if row_id not in self.happy_data:
                self.happy_data[row_id] = {"sample_id": samp}
            for fn in rdr.fieldnames:
                self.happy_data[row_id][fn] = row[fn]
                if 'METRIC' in fn:
                    gs_id = "{}/{}_{}".format(row["Type"], row["Filter"], fn)
                    self.general_stats_data[samp][gs_id] = row[fn]

        return 1


def gen_headers():
    h = OrderedDict()
    h["METRIC.Recall"] = {
        "title": "Recall",
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
    }

    h["TRUTH.TP"] = {
        "title": "Truth: True Positive",
        "description": "Number of true-positive calls in truth representation (counted via the truth sample column)",
        "format": None,
    }

    h["TRUTH.FN"] = {
        "title": "Truth: False Negative",
        "description": "Number of false-negative calls = calls in truth without matching query call",
        "format": None,
    }

    h["QUERY.TOTAL"] = {
        "title": "Query: Total",
        "description": "Total number of query calls",
        "format": None,
    }

    h["QUERY.TP"] = {
        "title": "Query: True Positive",
        "description": "Number of true positive calls in query representation (counted via the query sample column)",
        "format": None,
    }

    h["QUERY.FP"] = {
        "title": "Query: False Positive",
        "description": "Number of false-positive calls in the query file (mismatched query calls within the confident regions)",
        "format": None,
    }

    h["QUERY.UNK"] = {
        "title": "Query: Unknown",
        "description": "Number of query calls outside the confident regions",
        "format": None,
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

    return h
