""" MultiQC module to parse output from MultiVCFAnalyzer """


import json
import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """MultiVCFAnalyzer module"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="MultiVCFAnalyzer",
            anchor="multivcfanalyzer",
            href="https://github.com/alexherbig/MultiVCFAnalyzer",
            info="""combines multiple VCF files in a coherent way,
            can produce summary statistics and downstream analysis formats for phylogeny reconstruction.""",
            doi="10.1038/nature13591",
        )

        # Find and load any MultiVCFAnalyzer reports
        self.mvcf_data = dict()

        # Find and load JSON file
        for f in self.find_log_files("multivcfanalyzer", filehandles=True):
            self.parse_data(f)

        # Filter samples
        self.mvcf_data = self.ignore_samples(self.mvcf_data)

        # Return if no samples found
        if len(self.mvcf_data) == 0:
            raise ModuleNoSamplesFound

        # Add in extra columns to data file
        self.compute_perc_hets()

        # Add in extra stuff like Snp(hom) to data
        self.computeSnpHom()

        # Save data output file
        self.write_data_file(self.mvcf_data, "multiqc_multivcfanalyzer")

        # Add to General Statistics
        self.addSummaryMetrics()

        # Add MultiVCFAnalyzer Table Section
        self.add_section(name="Summary metrics", anchor="mvcf_table", plot=self.addTable())

        # Add MultiVCFAnalyzer Barplot Section
        self.add_section(
            name="Call statistics barplot",
            anchor="mvcf-barplot",
            helptext="""
        MultiVCFAnalyzer has a filtering step during the merge, where SNPs of low quality are discarded.
        This plot shows the number of SNPs that fell in to the different MultiVCFAnalyzer categories:

        * _SNP Calls (all):_ Total number of non-reference homzygous and heterozygous calls made
        * _Filtered SNP call:_ Number of non-reference positions excluded by user-supplied list.
        * _Number of Reference Calls:_ Number of reference calls made
        * _Discarded Reference Calls:_ Number of reference positions not reaching genotyping or coverage thresholds
        * _Discarded SNP Call:_ Number of non-reference positions not reaching enough coverage.
        * _No Call:_ Number of positions with no call made as reported by GATK
        * _Unhandled Genotypes:_ Number of positions where more than two possible alleles occured and were discarded
        """,
            plot=self.addBarplot(),
        )

    def parse_data(self, f):
        try:
            data = json.load(f["f"])
        except Exception as e:
            log.debug(e)
            log.warning(f"Could not parse MultiVCFAnalyzer JSON: '{f['fn']}'")
            return

        version = data.get("metadata", {}).get("version", None)

        # Parse JSON data to a dict
        for s_name, metrics in data.get("metrics", {}).items():
            s_clean = self.clean_s_name(s_name, f)
            if s_clean in self.mvcf_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_clean}")

            if version is not None:
                self.add_software_version(version, s_clean)

            self.add_data_source(f, s_clean)
            self.mvcf_data[s_clean] = dict()
            for snp_prop, value in metrics.items():
                self.mvcf_data[s_clean][snp_prop] = value

    # Compute % heterozygous snp alleles and add to data
    def compute_perc_hets(self):
        """Take the parsed stats from MultiVCFAnalyzer and add one column per sample"""
        for sample in self.mvcf_data:
            try:
                self.mvcf_data[sample]["Heterozygous SNP alleles (percent)"] = (
                    self.mvcf_data[sample]["SNP Calls (het)"] / self.mvcf_data[sample]["SNP Calls (all)"]
                ) * 100
            except ZeroDivisionError:
                self.mvcf_data[sample]["Heterozygous SNP alleles (percent)"] = "NA"

    def computeSnpHom(self):
        """Computes snp(hom) for data present by MultiVCFAnalyzer"""
        for sample in self.mvcf_data:
            self.mvcf_data[sample]["SNP Calls (hom)"] = (
                (self.mvcf_data[sample]["SNP Calls (all)"]) - self.mvcf_data[sample]["SNP Calls (het)"]
            )

    def addSummaryMetrics(self):
        """Take the parsed stats from MultiVCFAnalyzer and add it to the main plot"""

        headers = {
            "SNP Calls (all)": {
                "title": "SNPs",
                "description": "Total number of non-reference calls",
                "scale": "RdBu",
                "shared_key": "snp_call",
                "format": "{:,.0f}",
            },
            "SNP Calls (hom)": {
                "title": "Hom SNPs",
                "description": "Total number of non-reference calls passing homozygosity thresholds",
                "scale": "RdYlGn",
                "hidden": True,
                "shared_key": "snp_call",
                "format": "{:,.0f}",
            },
            "SNP Calls (het)": {
                "title": "Het SNPs",
                "description": "Total number of non-reference calls not passing homozygosity thresholds",
                "scale": "YlGn",
                "hidden": True,
                "shared_key": "snp_call",
                "format": "{:,.0f}",
            },
            "Heterozygous SNP alleles (percent)": {
                "title": "% Hets",
                "description": "Percentage of heterozygous SNP alleles",
                "scale": "OrRd",
                "max": 100,
                "min": 0,
            },
        }

        self.general_stats_addcols(self.mvcf_data, headers)

    def addTable(self):
        """Take the parsed stats from MultiVCFAnalyzer and add it to the MVCF Table"""
        headers = {
            "allPos": {
                "title": "Bases in Final Alignment",
                "description": "Length of FASTA file in base pairs (bp)",
                "scale": "BrBG",
                "shared_key": "calls",
                "format": "{:,.0f}",
            },
            "SNP Calls (all)": {
                "title": "SNPs",
                "description": "Total number of non-reference calls made",
                "scale": "OrRd",
                "shared_key": "snp_call",
                "format": "{:,.0f}",
            },
            "Heterozygous SNP alleles (percent)": {
                "title": "% Hets",
                "description": "Percentage of heterozygous SNP alleles",
                "scale": "PuBu",
                "max": 100,
                "min": 0,
            },
            "SNP Calls (hom)": {
                "title": "Hom SNPs",
                "description": "Total number of non-reference calls passing homozygosity thresholds",
                "scale": "RdYlGn",
                "shared_key": "snp_call",
                "format": "{:,.0f}",
            },
            "SNP Calls (het)": {
                "title": "Het SNPs",
                "description": "Total number of non-reference calls not passing homozygosity thresholds",
                "scale": "RdYlGn",
                "shared_key": "snp_call",
                "format": "{:,.0f}",
            },
            "discardedVarCall": {
                "title": "Discarded SNP Call",
                "description": "Number of non-reference positions not reaching genotyping or coverage thresholds",
                "scale": "PuRd",
                "shared_key": "calls",
                "format": "{:,.0f}",
            },
            "filteredVarCall": {
                "title": "Filtered SNP Call",
                "description": "Number of positions ignored defined in user-supplied filter list",
                "scale": "RdGy",
                "shared_key": "calls",
                "format": "{:,.0f}",
            },
            "refCall": {
                "title": "Reference Calls",
                "description": "Number of reference calls made",
                "scale": "Spectral",
                "shared_key": "calls",
                "format": "{:,.0f}",
            },
            "discardedRefCall": {
                "title": "Discarded Reference Call",
                "description": "Number of reference positions not reaching genotyping or coverage thresholds",
                "scale": "YlGnBu",
                "shared_key": "calls",
                "format": "{:,.0f}",
            },
            "coverage (fold)": {
                "title": "Average Call Coverage",
                "description": "Average number of reads covering final calls",
                "scale": "OrRd",
                "shared_key": "coverage",
                "suffix": "X",
            },
            "coverage (percent)": {
                "title": "% Reference with Calls",
                "description": "Percent coverage of all positions with final calls",
                "scale": "PuBuGn",
                "shared_key": "coverage",
                "suffix": "%",
                "max": 100,
                "min": 0,
            },
            "unhandledGenotype": {
                "title": "Unhandled Genotypes",
                "description": "Number of positions discarded due to presence of more than one alternate allele",
                "scale": "BuPu",
                "shared_key": "snp_count",
                "format": "{:,.0f}",
            },
            "noCall": {
                "title": "Positions with No Call",
                "description": "Number of positions with no call made as reported by GATK",
                "scale": "GnBu",
                "shared_key": "calls",
                "format": "{:,.0f}",
            },
        }

        # Separate table config
        table_config = {
            "namespace": "MultiVCFAnalyzer",  # Name for grouping. Prepends desc and is in Config Columns modal
            "id": "mvcf-table",  # ID used for the table
            "table_title": "MultiVCFAnalyzer Results",  # Title of the table. Used in the column config modal
        }
        tab = table.plot(self.mvcf_data, headers, table_config)
        return tab

    def addBarplot(self):
        """Take the parsed stats from MultiVCFAnalyzer and add it to the MVCF Table"""
        # SNP Calls (all): Green CHECK
        # Number of Reference Calls: Blue CHECK
        # Discarded SNP Call: Orange CHECK
        # Discarded Reference Call: Red [if avaliable] CHECK
        # Positions with No Call: Black CHECK
        cats = {
            "SNP Calls (hom)": {"name": "SNP Calls (hom)", "color": "#01665e"},
            "SNP Calls (het)": {"name": "SNP Calls (het)", "color": "#5ab4ac"},
            "refCall": {"name": "Reference Calls", "color": "#9ebcda"},
            "discardedVarCall": {"name": "Discarded SNP Call", "color": "#f7a35c"},
            "filteredVarCall": {
                "name": "Filtered SNP Call",
            },
            "discardedRefCall": {"name": "Discarded Reference Call", "color": "#e34a33"},
            "unhandledGenotype": {"name": "Positions with an unknown Genotype Call", "color": "#252525"},
            "noCall": {"name": "Positions with No Call", "color": "#252525"},
        }

        config = {
            # Building the plot
            "id": "mvcf_barplot",  # HTML ID used for plot
            "hide_zero_cats": True,  # Hide categories where data for all samples is 0
            # Customising the plot
            "title": "MultiVCFAnalyzer: Call Categories",  # Plot title - should be in format "Module Name: Plot Title"
            "ylab": "Total # Positions",  # X axis label
            "xlab": None,
            "stacking": "normal",  # Set to None to have category bars side by side
            "use_legend": True,  # Show / hide the legend
        }
        return bargraph.plot(self.mvcf_data, cats, config)
