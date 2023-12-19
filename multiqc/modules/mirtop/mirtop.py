""" MultiQC module to parse output from mirtop"""


import json
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="mirtop",
            anchor="mirtop",
            href="https://github.com/miRTop/mirtop/",
            info="is a command line tool to annotate miRNAs and isomiRs and compute general statistics using the mirGFF3 format.",
            doi="10.5281/zenodo.45385",  # Zenodo won't load this page for me as I write this, but it's the listed DOI.
        )

        # Find and load any mirtop reports
        self.mirtop_data = dict()
        for f in self.find_log_files("mirtop"):
            self.parse_mirtop_report(f)
            self.add_data_source(f)

        # Filter out ignored samples (given with --ignore-samples option)
        self.mirtop_data = self.ignore_samples(self.mirtop_data)

        # Raise error if dict is empty
        if len(self.mirtop_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.mirtop_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.mirtop_data, "multiqc_mirtop")

        # Helper variables
        self.general_helptext = """
            The different isomiR types are:

            * `iso_3p`: a sequence with a 3' end difference because of trimming or templated tailing
            * `iso_5p`: a sequence with a 5' end difference because of trimming or templated tailing
            * `iso_add3p`: a sequence with non templated tailing in the 3' end
            * `iso_add5p`: a sequence with non templated tailing in the 5' end
            * `iso_snv`: a sequence with a single nucleotide variant

            The `ref_miRNA` label corresponds to the reference miRNA (canonical sequence).

        """
        self.isomir_cats = [
            "ref_miRNA",
            "iso_3p",
            "iso_5p",
            "iso_add3p",
            "iso_add5p",
            "iso_add",
            "iso_snv",
            "iso_snv_seed",
            "iso_snv_central_offset",
            "iso_snv_central",
            "iso_snv_central_supp",
            "iso_snp",
            "iso_snp_seed",
            "iso_snp_central_offset",
            "iso_snp_central",
            "iso_snp_central_supp",
        ]

        # Calculate aggregate iso_snp counts
        self.mirtop_data_snp_aggregate = self.aggregate_snps_in_samples()

        # Create summary table
        self.mirtop_stats_table()

        # Create detailed plots
        self.mirtop_read_count()
        self.mirtop_unique_read_count()
        self.mirtop_mean_read_count()

    def parse_mirtop_report(self, f):
        """Parse the mirtop log file."""

        content = json.loads(f["f"])
        version = content.get("meta", {}).get("version")

        for s_name in content.get("metrics", {}).keys():
            cleaned_s_name = self.clean_s_name(s_name, f)
            ## Check for sample name duplicates
            if cleaned_s_name in self.mirtop_data:
                log.debug(f"Duplicate sample name found! Overwriting: {cleaned_s_name}")
            parsed_data = content["metrics"][s_name]
            # Sum the isomiR and ref_miRNA counts if present.
            parsed_data["read_count"] = parsed_data.get("isomiR_sum", 0) + parsed_data.get("ref_miRNA_sum", 0)
            if parsed_data["read_count"] > 0:
                parsed_data["isomiR_perc"] = (parsed_data.get("isomiR_sum", 0) / parsed_data.get("read_count", 0)) * 100
            else:
                parsed_data["isomiR_perc"] = 0.0
            self.mirtop_data[cleaned_s_name] = parsed_data

            if version is not None:
                version = version.strip("v")
                self.add_software_version(version, cleaned_s_name)

    def aggregate_snps_in_samples(self):
        """Aggregate info for iso_snp isomiRs (for clarity). "Mean" section will be recomputed"""
        snv_aggr = {}  ## sub dict with all infos except for snps
        for sample in self.mirtop_data:
            snv_aggr[sample] = {
                key: self.mirtop_data[sample][key] for key in self.mirtop_data[sample] if "iso_snp" not in key
            }
            snv_aggr[sample]["iso_snv_sum"] = sum(
                [self.mirtop_data[sample][key] for key in self.mirtop_data[sample] if "iso_snp" in key and "sum" in key]
            )
            snv_aggr[sample]["iso_snv_count"] = sum(
                [
                    self.mirtop_data[sample][key]
                    for key in self.mirtop_data[sample]
                    if "iso_snp" in key and "count" in key
                ]
            )

            if snv_aggr[sample]["iso_snv_count"] > 0:
                snv_aggr[sample]["iso_snv_mean"] = snv_aggr[sample]["iso_snv_sum"] / snv_aggr[sample]["iso_snv_count"]
            else:
                snv_aggr[sample]["iso_snv_mean"] = 0
        return snv_aggr

    def mirtop_stats_table(self):
        """Take the parsed stats from the mirtop report and add them to the
        basic stats table at the top of the report"""

        headers = {
            "ref_miRNA_sum": {
                "title": f"{config.read_count_prefix} Ref miRNA reads",
                "description": f"Read counts summed over all reference miRNAs ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "scale": "PuBu",
            },
            "isomiR_perc": {
                "title": "IsomiR %",
                "description": "% of total read counts corresponding to isomiRs",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "YlOrRd",
            },
            "isomiR_sum": {
                "title": f"{config.read_count_prefix} IsomiR reads",
                "description": f"Read counts summed over all isomiRs ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "scale": "Oranges",
            },
            "read_count": {
                "title": f"{config.read_count_prefix} Reads",
                "description": "Total read counts - both isomiRs and reference miRNA ({})".format(
                    config.read_count_desc
                ),
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "scale": "BuGn",
            },
        }

        self.general_stats_addcols(self.mirtop_data, headers)

    def filter_plot_data(self, plot_key):
        """Helper function to return the subset of data that has been requested for this plot"""
        plot_data = {}
        for s_name in self.mirtop_data_snp_aggregate:
            plot_data[s_name] = {}
            for key, data in self.mirtop_data_snp_aggregate[s_name].items():
                if plot_key in key and "isomiR" not in key and "read" not in key:
                    plot_data[s_name][key] = data
        return plot_data

    def get_plot_cats(self, plot_type):
        """Return the plot categories for the given plot"""
        cats_section = dict()
        for base_key in self.isomir_cats:
            cat_key = f"{base_key}_{plot_type}"
            cats_section[cat_key] = {"name": base_key}
        return cats_section

    def mirtop_read_count(self):
        """Generate barplot for the read count plot"""
        p_config = {"id": "mirtop_read_count_plot", "title": "mirtop: IsomiR read counts", "ylab": "Read counts"}
        self.add_section(
            name="IsomiR read counts",
            description="""
            Total counts of reads aligned for each isomiR type, over all detected miRNAs.
            """,
            helptext="""
            The total counts of reads detected as reference miRNA sequences is also shown.

            Since a read can belong to 2 (or more) different isomiRs types (e.g `iso_3p` and `iso_5p`),
            the cumulative read counts shown in this plot for a sample can be higher than its total
            read count shown in the general statistics.

            For each sample, the mean counts of each type of isomiRs over all detected
            miRNAs is displayed in a different color.
            """
            + self.general_helptext,
            plot=bargraph.plot(self.filter_plot_data("sum"), self.get_plot_cats("sum"), p_config),
        )

    def mirtop_unique_read_count(self):
        """Generate the section for the Unique Read Count plot"""
        p_config = {
            "id": "mirtop_unique_read_count_plot",
            "title": "mirtop: IsomiR unique read counts",
            "ylab": "Unique sequences",
        }
        self.add_section(
            name="IsomiR unique read counts",
            anchor="mirtop_unique_read_count",
            description="""
            The number of distinct sequences detected for each isomiR type, over all miRNAs.
            """,
            helptext="""
            The number of reference miRNA sequences detected is also shown.

            For each sample, the number of miRNAs with each type of isomiRs, is displayed in a different color.
            """
            + self.general_helptext,
            plot=bargraph.plot(self.filter_plot_data("count"), self.get_plot_cats("count"), p_config),
        )

    def mirtop_mean_read_count(self):
        """Generate the section for the Mean Read Count plot"""
        p_config = {"id": "mirtop_mean_read_count_plot", "title": "mirtop: Mean isomiR read counts", "ylab": "Means"}
        self.add_section(
            name="Mean isomiR read counts",
            anchor="mirtop_mean_read_count",
            description="""
            Mean counts for each isomiR type, over all detected miRNAs.
            """,
            helptext="""
            The mean counts of reads detected as reference miRNA sequences is also shown.

            For each sample, the mean counts of each type of isomiRs over all detected miRNAs is displayed in a different color.
            """
            + self.general_helptext,
            plot=bargraph.plot(self.filter_plot_data("mean"), self.get_plot_cats("mean"), p_config),
        )
