""" MultiQC module to parse output from SexdetErrmine """


import json
import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, scatter

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """SexDeterrmine module"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="SexDetErrmine",
            anchor="sexdeterrmine",
            href="https://github.com/TCLamnidis/Sex.DetERRmine",
            info="""A python script to calculate the relative coverage of X and Y chromosomes,
            and their associated error bars, from the depth of coverage at specified SNPs.""",
            doi="10.1038/s41467-018-07483-5",
        )

        # Find and load any DeDup reports
        self.sexdet_data = dict()

        # Find and load JSON file
        for f in self.find_log_files("sexdeterrmine", filehandles=True):
            self.parse_data(f)

        # Filter samples
        self.sexdet_data = self.ignore_samples(self.sexdet_data)

        # Return if no samples found
        if len(self.sexdet_data) == 0:
            raise ModuleNoSamplesFound

        # Save data output file
        self.write_data_file(self.sexdet_data, "multiqc_sexdeter_metrics")

        # Add to General Statistics
        self.addSummaryMetrics()

        # Plots
        self.snp_rate_scatterplot()
        self.read_count_barplot()
        self.snp_count_barplot()

    def parse_data(self, f):
        try:
            data = json.load(f["f"])
        except Exception as e:
            log.debug(e)
            log.warning(f"Could not parse SexDeterrmine JSON: '{f['fn']}'")
            return

        # Get the version
        version = str(data["Metadata"]["version"])

        # Parse JSON data to a dict
        for s_name in data:
            if s_name == "Metadata":
                continue

            s_clean = self.clean_s_name(s_name, f)
            if s_clean in self.sexdet_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_clean}")

            self.add_data_source(f, s_clean)
            self.add_software_version(version, s_clean)
            self.sexdet_data[s_clean] = dict()

            for k, v in data[s_name].items():
                try:
                    self.sexdet_data[s_clean][k] = float(v)
                except ValueError:
                    self.sexdet_data[s_clean][k] = v

    def addSummaryMetrics(self):
        """Take the parsed stats from SexDetErrmine and add it to the main plot"""

        headers = {
            "RateErrX": {
                "title": "Err Rate X",
                "description": "Rate of Error for Chr X",
                "scale": "OrRd",
                "hidden": True,
                "shared_key": "snp_err_rate",
            },
            "RateErrY": {
                "title": "Err Rate Y",
                "description": "Rate of Error for Chr Y",
                "scale": "OrRd",
                "hidden": True,
                "shared_key": "snp_err_rate",
            },
            "RateX": {
                "title": "Rate X",
                "description": "Number of positions on Chromosome X vs Autosomal positions.",
                "scale": "PuBuGn",
                "shared_key": "snp_count",
            },
            "RateY": {
                "title": "Rate Y",
                "description": "Number of positions on Chromosome Y vs Autosomal positions.",
                "scale": "BuPu",
                "shared_key": "snp_count",
            },
        }

        self.general_stats_addcols(self.sexdet_data, headers)

    def read_count_barplot(self):
        """Make a bar plot showing read counts on Autosomal, X and Y chr"""
        cats = {
            "NR Aut": {"name": "Autosomal Reads"},
            "NrX": {"name": "Reads on X"},
            "NrY": {"name": "Reads on Y"},
        }
        config = {
            "id": "sexdeterrmine-readcounts-plot",
            "title": "SexDetErrmine: Read Counts",
            "ylab": "# Reads",
        }

        self.add_section(
            name="Read Counts",
            anchor="sexdeterrmine-readcounts",
            description="The number of reads covering positions on the autosomes, X and Y chromosomes.",
            plot=bargraph.plot(self.sexdet_data, cats, config),
        )

    def snp_rate_scatterplot(self):
        """Make a scatter plot showing relative coverage on X and Y chr"""
        data = {}
        for sample in self.sexdet_data:
            try:
                data[sample] = {"x": self.sexdet_data[sample]["RateX"], "y": self.sexdet_data[sample]["RateY"]}
            except KeyError:
                pass

        config = {
            "id": "sexdeterrmine-rate-plot",
            "title": "SexDetErrmine: Relative coverage",
            "ylab": "Relative Cov. on Y",
            "xlab": "Relative Cov. on X",
        }

        if len(data) > 0:
            self.add_section(
                name="Relative Coverage",
                anchor="sexdeterrmine-rates",
                description="The coverage on the X vs Y chromosome, relative to coverage on the Autosomes.",
                helptext="""
                Males are expected to have a roughly equal X- and Y-rates, while females are expected to have a Y-rate of 0 and an X-rate of 1.
                Placement between the two clusters can be indicative of contamination, while placement with higher than expected X- and/or Y-rates can be indicative of sex chromosome aneuploidy.
                """,
                plot=scatter.plot(data, config),
            )

    def snp_count_barplot(self):
        """Make a bar plot showing read counts on Autosomal, X and Y chr"""
        cats = {
            "Snps Autosomal": {"name": "Autosomal SNPs"},
            "XSnps": {"name": "SNPs on X"},
            "YSnps": {"name": "SNPs on Y"},
        }

        config = {
            "id": "sexdeterrmine-snps-plot",
            "title": "SexDetErrmine: SNP Counts",
            "ylab": "# Reads",
        }

        self.add_section(
            name="SNP Counts",
            anchor="sexdeterrmine-snps",
            description="Total number of SNP positions. When supplied with a BED file, this includes only positions specified there.",
            plot=bargraph.plot(self.sexdet_data, cats, config),
        )
