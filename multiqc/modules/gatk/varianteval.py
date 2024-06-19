"""MultiQC submodule to parse output from GATK varianteval"""

import logging

from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class VariantEvalMixin:
    def parse_gatk_varianteval(self):
        """Find GATK varianteval logs and parse their data"""

        self.gatk_varianteval = dict()
        for f in self.find_log_files("gatk/varianteval", filehandles=True):
            parsed_data = parse_single_report(f["f"])
            if len(parsed_data) > 1:
                if f["s_name"] in self.gatk_varianteval:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                self.add_data_source(f, section="varianteval")
                self.gatk_varianteval[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.gatk_varianteval = self.ignore_samples(self.gatk_varianteval)

        n_reports_found = len(self.gatk_varianteval)
        if n_reports_found == 0:
            return 0

        log.info(f"Found {n_reports_found} VariantEval reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file (restructure first)
        self.write_data_file(self.gatk_varianteval, "multiqc_gatk_varianteval")

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
            "description": f"TiTV ratio from variants found in '{titv_ref}'",
            "min": 0,
            "scale": "Blues",
            "shared_key": "titv_ratio",
        }
        varianteval_headers["novel_titv"] = {
            "title": "TiTV ratio (novel)",
            "description": f"TiTV ratio from variants NOT found in '{titv_ref}'",
            "min": 0,
            "scale": "Blues",
            "shared_key": "titv_ratio",
        }
        varianteval_headers["called_titv"] = {
            "title": "TiTV ratio (called)",
            "description": f"TiTV ratio from variants found in '{titv_ref}'",
            "min": 0,
            "scale": "Blues",
            "shared_key": "titv_ratio",
        }
        varianteval_headers["filtered_titv"] = {
            "title": "TiTV ratio (filtered)",
            "description": f"TiTV ratio from variants NOT found in '{titv_ref}'",
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


def parse_single_report(f):
    """Parse a gatk varianteval varianteval"""
    # Fixme: Separate GATKReport parsing and data subsetting. A GATKReport parser now available from the GATK MultiqcModel.

    data = dict()
    in_CompOverlap = False
    in_CountVariants = False
    in_TiTv = False
    for line in f:
        # Detect section headers
        if "#:GATKTable:CompOverlap" in line:
            in_CompOverlap = True
        elif "#:GATKTable:CountVariants" in line:
            in_CountVariants = True
        elif "#:GATKTable:TiTvVariantEvaluator" in line:
            in_TiTv = True
        else:
            # Parse contents using nested loops
            if in_CompOverlap:
                headers = line.split()
                while in_CompOverlap:
                    line = f.readline().strip("\n")
                    if len(line) < len(headers):
                        in_CompOverlap = False
                        break
                    d = dict()
                    for i, s in enumerate(line.split()):
                        d[headers[i]] = s
                    if d.get("Novelty") == "all" or d.get("Filter") == "raw":
                        if "CompFeatureInput" in d:
                            data["reference"] = d["CompFeatureInput"]
                        if "CompRod" in d:
                            data["reference"] = d["CompRod"]
                        if "compRate" in d:
                            data["comp_rate"] = float(d["compRate"])
                        if "concordantRate" in d:
                            data["concordant_rate"] = float(d["concordantRate"])
                        if "nEvalVariants" in d:
                            data["eval_variants"] = int(d["nEvalVariants"])
                        if "novelSites" in d:
                            data["novel_sites"] = int(d["novelSites"])
                    elif d.get("Novelty") == "known":
                        if "nEvalVariants" in d:
                            data["known_sites"] = int(d["nEvalVariants"])

            elif in_CountVariants:
                headers = line.split()
                while in_CountVariants:
                    line = f.readline().strip("\n")
                    if len(line) < len(headers):
                        in_CountVariants = False
                        break
                    d = dict()
                    for i, s in enumerate(line.split()):
                        d[headers[i]] = s
                    if d.get("Novelty") == "all" or d.get("Filter") == "raw":
                        if "nSNPs" in d:
                            data["snps"] = int(d["nSNPs"])
                        if "nMNPs" in d:
                            data["mnps"] = int(d["nMNPs"])
                        if "nInsertions" in d:
                            data["insertions"] = int(d["nInsertions"])
                        if "nDeletions" in d:
                            data["deletions"] = int(d["nDeletions"])
                        if "nComplex" in d:
                            data["complex"] = int(d["nComplex"])
                        if "nSymbolic" in d:
                            data["symbolic"] = int(d["nSymbolic"])
                        if "nMixed" in d:
                            data["mixed"] = int(d["nMixed"])
                        if "nNoCalls" in d:
                            data["nocalls"] = int(d["nNoCalls"])
            elif in_TiTv:
                headers = line.split()
                data["titv_reference"] = "unknown"
                while in_TiTv:
                    line = f.readline().strip("\n")
                    if len(line) < len(headers):
                        in_TiTv = False
                        break
                    d = dict()
                    for i, s in enumerate(line.split()):
                        d[headers[i]] = s
                    if d.get("Novelty") == "known":
                        if "CompFeatureInput" in d:
                            data["titv_reference"] = d["CompFeatureInput"]
                        if "CompRod" in d:
                            data["titv_reference"] = d["CompRod"]
                        if "tiTvRatio" in d:
                            data["known_titv"] = float(d["tiTvRatio"])
                    elif d.get("Novelty") == "novel":
                        if "tiTvRatio" in d:
                            data["novel_titv"] = float(d["tiTvRatio"])
                    elif d.get("Filter") == "called":
                        if "CompFeatureInput" in d:
                            data["titv_reference"] = d["CompFeatureInput"]
                        if "CompRod" in d:
                            data["titv_reference"] = d["CompRod"]
                        if "tiTvRatio" in d:
                            data["called_titv"] = float(d["tiTvRatio"])
                    elif d.get("Filter") == "filtered":
                        if "CompFeatureInput" in d:
                            data["titv_reference"] = d["CompFeatureInput"]
                        if "CompRod" in d:
                            data["titv_reference"] = d["CompRod"]
                        if "tiTvRatio" in d:
                            data["filtered_titv"] = float(d["tiTvRatio"])

    return data


def count_variants_barplot(data):
    """Return HTML for the Variant Counts barplot"""
    keys = {
        "snps": {"name": "SNPs"},
        "mnps": {"name": "MNPs"},
        "insertions": {"name": "Insertions"},
        "deletions": {"name": "Deletions"},
        "complex": {"name": "Complex"},
        "symbolic": {"name": "Symbolic"},
        "mixed": {"name": "Mixed"},
        "nocalls": {"name": "No-calls"},
    }

    plot_conf = {
        "id": "gatk_varianteval_variant_plot",
        "title": "GATK VariantEval: Variant Counts",
        "ylab": "# Variants",
        "cpswitch_counts_label": "Number of Variants",
    }
    return bargraph.plot(data, keys, plot_conf)


def comp_overlap_table(data):
    """Build a table from the comp overlaps output."""
    headers = {
        "comp_rate": {
            "title": "Compare rate",
            "description": "Ratio of known variants found in the reference set.",
            "namespace": "GATK",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "Blues",
        },
        "concordant_rate": {
            "title": "Concordant rate",
            "description": "Ratio of variants matching alleles in the reference set.",
            "namespace": "GATK",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "Blues",
        },
        "eval_variants": {
            "title": "M Evaluated variants",
            "description": "Number of called variants (millions)",
            "namespace": "GATK",
            "min": 0,
            "modify": lambda x: float(x) / 1000000.0,
        },
        "known_sites": {
            "title": "M Known sites",
            "description": "Number of known variants (millions)",
            "namespace": "GATK",
            "min": 0,
            "modify": lambda x: float(x) / 1000000.0,
        },
        "novel_sites": {
            "title": "M Novel sites",
            "description": "Number of novel variants (millions)",
            "namespace": "GATK",
            "min": 0,
            "modify": lambda x: float(x) / 1000000.0,
        },
    }
    table_html = table.plot(data, headers, {"id": "gatk_compare_overlap", "title": "GATK - Compare Overlap"})
    return table_html
