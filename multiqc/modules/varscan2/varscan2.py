""" MultiQC module to parse output files from VarScan2 """


import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="VarScan2",
            anchor="varscan",
            href="http://dkoboldt.github.io/varscan/",
            info="variant detection in massively parallel sequencing data",
            doi=["10.1101/gr.129684.111", "10.1093/bioinformatics/btp373"],
        )

        # Find and load VarScan2 reports - there are 3 different ones, but all with identical content (differentiated by header)
        self.varscan2_data = dict()
        for f in self.find_log_files("varscan2/mpileup2snp", filehandles=True):
            parsed_data = self.parse_varscan(f)
            s_name = self.clean_s_name(parsed_data["sample_name"], f)
            # Drop existing sample_name from dict now
            del parsed_data["sample_name"]
            if parsed_data is not None and len(parsed_data) > 0:
                self.varscan2_data[s_name] = parsed_data
                self.add_data_source(f, s_name, section="mpileup2snp")

        for f in self.find_log_files("varscan2/mpileup2indel", filehandles=True):
            parsed_data = self.parse_varscan(f)
            s_name = self.clean_s_name(parsed_data["sample_name"], f)
            # Drop existing sample_name from dict now
            del parsed_data["sample_name"]
            if parsed_data is not None and len(parsed_data) > 0:
                self.varscan2_data[s_name] = parsed_data
                self.add_data_source(f, s_name, section="mpileup2indel")

        for f in self.find_log_files("varscan2/mpileup2cns", filehandles=True):
            parsed_data = self.parse_varscan(f)
            s_name = self.clean_s_name(parsed_data["sample_name"], f)
            # Drop existing sample_name from dict now
            del parsed_data["sample_name"]
            if parsed_data is not None and len(parsed_data) > 0:
                self.varscan2_data[s_name] = parsed_data
                self.add_data_source(f, s_name, section="mpileup2cns")

        # Filter to strip out ignored sample names
        self.varscan2_data = self.ignore_samples(self.varscan2_data)

        # Warning when no files are found
        if len(self.varscan2_data) == 0:
            raise ModuleNoSamplesFound

        # Write parsed data to a file
        self.write_data_file(self.varscan2_data, "multiqc_varscan2_summary")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Found reports or not?
        log.info(f"Found {len(self.varscan2_data)} reports")

        # Basic Stats Table
        self.varscan2_general_stats_table()

        # Basic barplot section
        self.varscan2_counts_barplot()

    # Varscan2 reports in SNP mode only snps, in indel only indel and in CNS all variants found
    # Total variants = SNPs + Indels
    def parse_varscan(self, f):
        """Parse a VarScan2 report"""
        parsed = dict()
        regexes = {
            "sample_name": r"(?:Reading input from )(\w+.+)",
            "min_coverage": r"(?:Min coverage:)\s(\d+)",
            "min_reads2": r"(?:Min reads2:)\s(\d+)",
            "min_var_freq": r"(?:Min var freq:)\s(\d+\.\d+)",
            "min_avg_qual": r"(?:Min avg qual:)\s(\d+)",
            "p_val_threshold": r"(?:P-value thresh:)\s(\d+\.\d+)",
            "bases_in_pileup": r"(\d+)(?:\sbases in pileup file)",
            "variant_total": r"(\d+)(?:\svariant positions \()",
            "variant_snps": r"(?:\d+)(?:\svariant positions \()(\d+)",
            "variant_indels": r"(?:\d+)(?:\svariant positions \()(?:\d+ SNP, )(\d+)",
            "failed_strand_filter": r"(\d+)(?:\swere failed by the strand-filter)",
            "variant_reported_total": r"(\d+)(?:\svariant positions reported)",
            "variant_reported_snps": r"(?:\d+)(?:\svariant positions reported \()(\d+)",
            "variant_reported_indels": r"(?:\d+)(?:\svariant positions reported \()(?:\d+ SNP, )(\d+)",
        }
        for line in f["f"]:
            # Search regexes for stats
            for k, r in regexes.items():
                match = re.search(r, line)
                if match:
                    if k not in ["sample_name", "p_val_threshold", "min_var_freq"]:
                        parsed[k] = int(match.group(1))
                    if k == "sample_name":
                        parsed[k] = match.group(1)
                    if k in ["p_val_threshold", "min_var_freq"]:
                        parsed[k] = float(match.group(1))

        # The failed_strand_filter combines the different classes, so recompute
        if "variant_total" in parsed and "variant_reported_total" in parsed:
            parsed["variant_total_failed"] = parsed["variant_total"] - parsed["variant_reported_total"]
        if "variant_snps" in parsed and "variant_reported_snps" in parsed:
            parsed["variant_snps_failed"] = parsed["variant_snps"] - parsed["variant_reported_snps"]
        if "variant_indels" in parsed and "variant_reported_indels" in parsed:
            parsed["variant_indels_failed"] = parsed["variant_indels"] - parsed["variant_reported_indels"]

        return parsed

    # Add to general stats table
    def varscan2_general_stats_table(self):
        """Take the parsed stats from the VarScan2 report and add it to the
        General Statistics table
        """

        headers = {
            "variant_reported_snps": {
                "title": "SNPs reported",
                "description": "Total number of SNPs reported.",
                "min": 0,
                "scale": "Spectral",
                "format": "{:,.0f}",
                "shared_key": "snps",
            },
            "variant_reported_indels": {
                "title": "INDELs reported",
                "description": "Total number INDELs reported.",
                "min": 0,
                "scale": "Spectral",
                "format": "{:,.0f}",
                "shared_key": "indels",
            },
            "variant_reported_total": {
                "title": "Variants reported",
                "description": "Total number of variants reported.",
                "min": 0,
                "scale": "Spectral",
                "format": "{:,.0f}",
                "shared_key": "variants",
                "hidden": True,
            },
            "variant_snps": {
                "title": "SNPs",
                "description": "Total number of SNPs detected",
                "min": 0,
                "scale": "BrBG",
                "format": "{:,.0f}",
                "shared_key": "snps",
                "hidden": True,
            },
            "variant_indels": {
                "title": "INDELs",
                "description": "Total number of INDELs detected",
                "min": 0,
                "scale": "BrBG",
                "format": "{:,.0f}",
                "shared_key": "indels",
                "hidden": True,
            },
            "variant_total": {
                "title": "Variants",
                "description": "Total number of variants detected",
                "min": 0,
                "scale": "BrBG",
                "format": "{:,.0f}",
                "shared_key": "variants",
                "hidden": True,
            },
            "failed_strand_filter": {
                "title": "Failed Strand Filter",
                "description": "Total number variants failing the strand-filter.",
                "min": 0,
                "scale": "YlOrBr",
                "format": "{:,.0f}",
                "shared_key": "variants",
                "hidden": True,
            },
            "bases_in_pileup": {
                "title": f"{config.base_count_prefix} Bases in Pileup",
                "description": "Number of bases in pileup input for VarScan2 ()",
                "scale": "Greens",
                "modify": lambda x: x * config.base_count_multiplier,
                "shared_key": "base_count",
                "hidden": True,
            },
        }

        self.general_stats_addcols(self.varscan2_data, headers)

    def varscan2_counts_barplot(self):
        # 146 variant positions (106 SNP, 40 indel)
        # 12 were failed by the strand-filter
        # 99 variant positions reported (99 SNP, 0 indel)

        # Specify the order of the different possible categories
        cats = [
            {
                "variant_reported_snps": {"name": "SNPs reported", "color": "#7fc97f"},
                "variant_snps_failed": {"name": "SNPs not reported", "color": "#fb8072"},
            },
            {
                "variant_reported_indels": {"name": "INDELs reported", "color": "#7fc97f"},
                "variant_indels_failed": {"name": "INDELs not reported", "color": "#fb8072"},
            },
        ]
        # Config for the plot
        pconfig = {
            "id": "varscan2_variant_counts_plot",
            "title": "VarScan2: Variants detected",
            "ylab": "Number of SNPs",
            "cpswitch_counts_label": "Number of Variants",
            "hide_zero_cats": False,
            "data_labels": [{"name": "SNPs", "ylab": "Number of SNPs"}, {"name": "INDELs", "ylab": "Number of INDELs"}],
        }

        self.add_section(
            name="Variants detected",
            anchor="varscan2_variant_counts",
            description="""
            This plot shows the total number of variant positions, broken down by those that were reported or not.
            """,
            plot=bargraph.plot([self.varscan2_data, self.varscan2_data], cats, pconfig),
        )
