# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from GATK tool AnalyzeSaturationMutagenesis """

import logging
from collections import OrderedDict

from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class AnalyzeSaturationMutagenesisMixin:
    def parse_gatk_analyze_saturation_mutagenesis(self):
        """Find GATK AnalyzeSaturationMutagenesis logs and parse their data"""

        self.gatk_analyze_saturation_mutagenesis = dict()
        for f in self.find_log_files("gatk/analyze_saturation_mutagenesis", filehandles=True):
            


            parsed_data = parse_readCount_file(f["f"])
            if len(parsed_data) > 1:
                if f["s_name"] in self.gatk_analyze_saturation_mutagenesis:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
                self.add_data_source(f, section="analyze_saturation_mutagenesis")
                self.gatk_analyze_saturation_mutagenesis[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.gatk_analyze_saturation_mutagenesis = self.ignore_samples(self.gatk_analyze_saturation_mutagenesis)

        n_reports_found = len(self.gatk_analyze_saturation_mutagenesis)
        if n_reports_found > 0:
            log.info("Found {} AnalyzeSaturationMutagenesis reports".format(n_reports_found))

            # Write parsed report data to a file (restructure first)
            self.write_data_file(self.gatk_analyze_saturation_mutagenesis, "multiqc_gatk_analyze_saturation_mutagenesis")

            # Get consensus TiTv references
            titv_ref = None
            for s_name in self.gatk_varianteval:
                if titv_ref is None:
                    titv_ref = self.gatk_varianteval[s_name]["titv_reference"]
                elif titv_ref != self.gatk_varianteval[s_name]["titv_reference"]:
                    titv_ref = "Multiple"
                    break

            # General Stats Table
            varianteval_headers = dict()
            varianteval_headers["known_titv"] = {
                "title": "TiTV ratio (known)",
                "description": "TiTV ratio from variants found in '{}'".format(titv_ref),
                "min": 0,
                "scale": "Blues",
                "shared_key": "titv_ratio",
            }
            varianteval_headers["novel_titv"] = {
                "title": "TiTV ratio (novel)",
                "description": "TiTV ratio from variants NOT found in '{}'".format(titv_ref),
                "min": 0,
                "scale": "Blues",
                "shared_key": "titv_ratio",
            }
            varianteval_headers["called_titv"] = {
                "title": "TiTV ratio (called)",
                "description": "TiTV ratio from variants found in '{}'".format(titv_ref),
                "min": 0,
                "scale": "Blues",
                "shared_key": "titv_ratio",
            }
            varianteval_headers["filtered_titv"] = {
                "title": "TiTV ratio (filtered)",
                "description": "TiTV ratio from variants NOT found in '{}'".format(titv_ref),
                "min": 0,
                "scale": "Blues",
                "shared_key": "titv_ratio",
            }
            self.general_stats_addcols(self.gatk_varianteval, varianteval_headers, "GATK VariantEval")

            # Variant Counts plot
            self.add_section(
                name="Variant Counts", anchor="gatk-count-variants", plot=count_variants_barplot(self.gatk_varianteval)
            )

            # Compare Overlap Table
            self.add_section(
                name="Compare Overlap", anchor="gatk-compare-overlap", plot=comp_overlap_table(self.gatk_varianteval)
            )

        return n_reports_found


def parse_readCounts_file(f):
    """Parse a readCounts output file from GATK AnalyzeSaturationMutagenesis
    These files are tab delimited, with some hierarchical structuring. """

    data = dict()

    in_disjoint = False
    in_overlapping = False

    for l in f:
        # Detect the two hierarchical sections
        if "disjoint pairs evaluated separately" in l:
            in_disjoint = True
            data['disjoint_total'], data['disjoint_total_percent'] = l.split('\t')[1:3]

        elif "overlapping pairs evaluated together" in l:
            in_overlapping = True
            data['overlapping_total'], data['overlapping_total_percent'] = l.split('\t')[1:3]

        # Parse the rest of the lines
        # The percentage is written as a string with a % character, so we parse it by stripping the '%' and casting to float.

        else:
            fields = l.split('\t')
            if in_disjoint:
                if fields[0] == '>>>Wild type:':
                    data['disjoint_wt'], data['disjoint_wt_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>>>Called variants:':
                    data['disjoint_called'], data['disjoint_called_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>>>Mate ignored:':
                    data['disjoint_mate_ignored'], data['disjoint_mate_ignored_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>>>Low quality variation:':
                    data['disjoint_low_q'], data['disjoint_low_q_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>>>Insufficient flank:':
                    data['disjoint_flank'], data['disjoint_flank_percent'] = fields[1], float(fields[2][:-1])
                    # Last line of disjoint pairs section, so unset flag
                    in_disjoint = False

            elif in_overlapping:
                if fields[0] == '>>>Wild type:':
                    data['overlapping_wt'], data['overlapping_wt_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>>>Called variants:':
                    data['overlapping_called'], data['overlapping_called_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>>>Mate ignored:':
                    data['overlapping_mate_ignored'], data['overlapping_mate_ignored_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>>>Low quality variation:':
                    data['overlapping_low_q'], data['overlapping_low_q_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>>>Insufficient flank:':
                    data['overlapping_flank'], data['overlapping_flank_percent'] = fields[1], float(fields[2][:-1])
                    # Last line of overlapping pairs section, so unset flag
                    in_overlapping = False
            
            else:
                if fields[0] == 'Total Reads:':
                    data['total'], data['total_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>Unmapped Reads:':
                    data['unmapped'], data['unmapped_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>LowQ Reads:':
                    data['lowq'], data['lowq_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>Evaluable Reads:':
                    data['evaluable'], data['evaluable_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == 'Total base calls:':
                    data['total_base_calls'], data['total_base_calls_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>Base calls evaluated for variants:':
                    data['evaluated_base_calls'], data['evaluated_base_calls_percent'] = fields[1], float(fields[2][:-1])
                elif fields[0] == '>Base calls unevaluated:':
                    data['unevaluated_base_calls'], data['unevaluated_base_calls_percent'] = fields[1], float(fields[2][:-1])

    return data

def count_variants_barplot(data):
    """Return HTML for the Variant Counts barplot"""
    keys = OrderedDict()
    keys["snps"] = {"name": "SNPs"}
    keys["mnps"] = {"name": "MNPs"}
    keys["insertions"] = {"name": "Insertions"}
    keys["deletions"] = {"name": "Deletions"}
    keys["complex"] = {"name": "Complex"}
    keys["symbolic"] = {"name": "Symbolic"}
    keys["mixed"] = {"name": "Mixed"}
    keys["nocalls"] = {"name": "No-calls"}

    plot_conf = {
        "id": "gatk_varianteval_variant_plot",
        "title": "GATK VariantEval: Variant Counts",
        "ylab": "# Variants",
        "cpswitch_counts_label": "Number of Variants",
    }
    return bargraph.plot(data, keys, plot_conf)


def comp_overlap_table(data):
    """Build a table from the comp overlaps output."""
    headers = OrderedDict()
    headers["comp_rate"] = {
        "title": "Compare rate",
        "description": "Ratio of known variants found in the reference set.",
        "namespace": "GATK",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "Blues",
    }
    headers["concordant_rate"] = {
        "title": "Concordant rate",
        "description": "Ratio of variants matching alleles in the reference set.",
        "namespace": "GATK",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "Blues",
    }
    headers["eval_variants"] = {
        "title": "M Evaluated variants",
        "description": "Number of called variants (millions)",
        "namespace": "GATK",
        "min": 0,
        "modify": lambda x: float(x) / 1000000.0,
    }
    headers["known_sites"] = {
        "title": "M Known sites",
        "description": "Number of known variants (millions)",
        "namespace": "GATK",
        "min": 0,
        "modify": lambda x: float(x) / 1000000.0,
    }
    headers["novel_sites"] = {
        "title": "M Novel sites",
        "description": "Number of novel variants (millions)",
        "namespace": "GATK",
        "min": 0,
        "modify": lambda x: float(x) / 1000000.0,
    }
    table_html = table.plot(data, headers, {"id": "gatk_compare_overlap", "table_title": "GATK - Compare Overlap"})
    return table_html
