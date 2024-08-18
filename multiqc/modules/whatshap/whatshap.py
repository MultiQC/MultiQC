import logging
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


ALL_CHROM = "ALL"


class MultiqcModule(BaseMultiqcModule):
    """
    The module is currently restricted to the output from `whatshap stats --tsv`.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="WhatsHap",
            anchor="whatshap",
            href="https://whatshap.readthedocs.io/",
            info="Phasing genomic variants using DNA reads (aka read-based phasing, or haplotype assembly)",
            doi="10.1101/085050",
        )

        # Store the whatshap stats results
        self.whatshap_stats = dict()

        # Iterate over all files we found
        for f in self.find_log_files("whatshap/stats", filehandles=True):
            sample, data = self.parse_whatshap_stats(f)
            self.whatshap_stats[sample] = data
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.whatshap_stats = self.ignore_samples(self.whatshap_stats)

        if not self.whatshap_stats:
            raise ModuleNoSamplesFound

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.whatshap_stats, "multiqc_whatshap_stats")

        # Add whatshap stats to general statistics
        self.whatshap_add_general_stats()

        # Create bargraph with total phased bp
        self.add_bargraph_total_phased()

        # Add a table with more extensive whatshap stats not shown in the
        # general statistics table
        self.add_stats_table()

    def parse_whatshap_stats(self, logfile):
        """Parse WhatsHap stats file"""

        def parse_numeric_values(data):
            """
            Process the data to convert it to the proper format

            In practice, this only involves the following changes:
                - Convert all numeric values to integers or floats
                - Replace 'nan' values with zero
            """
            for key, value in data.items():
                # Try to make the value an int
                try:
                    data[key] = int(data[key])
                except ValueError:
                    try:
                        data[key] = float(data[key])
                    except ValueError:
                        pass

                # Replace 'nan' with zero
                if value == "nan":
                    data[key] = None

        file_content = logfile["f"]
        # This will later be replaced by the sample name from the file, unless
        # we are parsing an empty file
        sample = logfile["s_name"]
        # Get the header and remove the "#" character from the first field
        header = next(file_content).strip().split("\t")
        header[0] = header[0].lstrip("#")

        # The results from this log file, a dictionary for every chromosome
        results = defaultdict(dict)
        # Parse the lines that have the data
        for line in file_content:
            fields = line.strip().split("\t")
            if len(header) != len(fields):
                continue

            data = {field: value for field, value in zip(header, fields)}
            if "chromosome" not in data or "sample" not in data:
                continue

            # Remove the sample and chromosome from the data
            sample = str(data.pop("sample"))
            chromosome = data.pop("chromosome")

            # Process the remaining data fields
            parse_numeric_values(data)

            # Insert the current line under chromosome
            results[chromosome] = data

        # If we were parsing an empty file, we return dummy data (basically all
        # zero) so the sample statistics are consistent, which eases parsing
        # later on. WhatsHap uses the fictional 'ALL' chromosome to store the
        # summary for all results, so we re-use that to store all zero's
        if not results:
            results = dict()
            results[ALL_CHROM] = defaultdict(int)

        # If ALL chromosome is not present, and there are multiple contigs available,
        # e.g. in case of a truncated input file, we summarize available contigs
        # as a fictional 'ALL' chromosome.
        if ALL_CHROM not in results and len(results) > 1:
            log.warning(
                f"Could not find the 'ALL' chromosome for {sample}. Make sure that the input is not truncated. MultiQC will summarize information for all available chromosomes instead."
            )
            sum_fields = {
                "variants",
                "phased",
                "unphased",
                "singletons",
                "blocks",
                "bp_per_block_sum",
                "heterozygous_variants",
                "heterozygous_snvs",
                "phased_snvs",
            }
            avg_fields = {
                "bp_per_block_avg",
                "variant_per_block_avg",
                "bp_per_block_avg",
            }
            results[ALL_CHROM] = defaultdict(int)
            for chrom in results:
                if chrom == ALL_CHROM:
                    continue
                for field in results[chrom]:
                    val = results[chrom][field]
                    if val:
                        if field in sum_fields:
                            results[ALL_CHROM][field] += int(results[chrom][field])
                        elif field in avg_fields:
                            results[ALL_CHROM][field] += float(results[chrom][field]) * results[chrom]["blocks"]
            for field in avg_fields:
                results[ALL_CHROM][field] /= results[ALL_CHROM]["blocks"]

        # Calculate the fraction of heterozygous variants that were phased
        for chrom, data in results.items():
            try:
                frac_het_phased = data["phased"] / data["heterozygous_variants"]
            except ZeroDivisionError:
                frac_het_phased = 0
            results[chrom]["frac_het_phased"] = frac_het_phased

        # Clean the sample name
        sample = self.clean_s_name(sample, logfile["root"])

        return sample, results

    @staticmethod
    def get_summary_field(sample_stats):
        """
        Return the summary field for sample_stats

        If there is only a single chromosome, we use that as the summary field.
        If there are multiple chromosomes, whatshap adds the 'ALL' field which
        contains the summary
        """
        # If there is only a single chromosome we use that as the summary
        # filed. Otherwise, we use the 'ALL' chromosome (it only gets added
        # by WhatsHap stats when there are multiple chromosomes).
        if len(sample_stats) == 1:
            summary_field = list(sample_stats)[0]
        else:
            summary_field = ALL_CHROM
        return summary_field

    def whatshap_add_general_stats(self):
        """Add WhatsHap stats to the general statistics table"""

        # Store the configuration for each column, this is also used to
        # determine which columns to add to the general statistics table
        # https://whatshap.readthedocs.io/en/latest/guide.html#the-tsv-statistics-format
        general_stats_headers = {
            "frac_het_phased": {
                "title": "% Phased Variants",
                "description": """Fraction of heterozygous variants
                                          that could be phased.
                                        """,
                "format": "{:,.0f}",
                "suffix": "%",
                "min": 0,
                "max": 100,
                "modify": lambda x: x * 100,
                "scale": "GnBu",
                "hidden": False,
            },
            "bp_per_block_avg": {
                "title": "Avg bp per Block",
                "description": """Description of the distribution of non-singleton
                                        block lengths, where the length of a block is the
                                        number of base pairs it covers minus 1. That is, a
                                        block with two variants at positions 2 and 5 has
                                        length 3.""",
                "format": "{:,.0f}",
                "scale": "Greens",
                "hidden": False,
            },
            "block_n50": {
                "title": "NG50",
                "description": """The NG50 value of the distribution of the block
                                        lengths. Interleaved blocks are cut in order to
                                        avoid artificially inflating this value.""",
                "format": "{:,.0f}",
                "scale": "Blues",
                "hidden": True,
            },
        }

        general = dict()
        for sample, sample_stats in self.whatshap_stats.items():
            # Get the summary field
            summary_field = self.get_summary_field(sample_stats)

            # The summary data we are adding to the general statistics
            summary_data = sample_stats[summary_field]

            general_stats = {field: summary_data[field] for field in general_stats_headers}

            general[sample] = general_stats

        self.general_stats_addcols(general, general_stats_headers)

    def add_bargraph_total_phased(self):
        """Add a bargraph of the total number of phased base pairs"""
        pdata = dict()
        for sample, sample_stats in self.whatshap_stats.items():
            # Get the summary field
            summary_field = self.get_summary_field(sample_stats)
            pdata[sample] = {"Phased Base Pairs": sample_stats[summary_field]["bp_per_block_sum"]}

        configuration = {
            "id": "multiqc_whatshap_phased_bp_plot",
            "title": "WhatsHap: Phased Basepairs per Sample",
            "ylab": "Base Pairs",
            "cpswitch": False,
        }

        keys = {"Phased Base Pairs": {"name": "Phased Base Pairs"}}

        # If the Phased Base Pairs is zero for all samples, we do not add the
        # WhatsHap section
        for sample, values in pdata.items():
            if values["Phased Base Pairs"]:
                # If we found one sample with data, we break out of the loop
                # and add the WhatsHap section
                break
        # If we completed the for loop, there is no data, so we return without
        # adding the WhatsHap section
        else:
            return

        self.add_section(
            name="Phased Basepairs per Sample",
            anchor="multiqc_whatshap_phased_bp",
            description="""
                    This plot show the total number of phased base pairs for
                    each sample.
                """,
            plot=bargraph.plot(pdata, keys, configuration),
        )

    def add_stats_table(self):
        """Add WhatsHap stats to the general statistics table"""

        # Store the configuration for each column, this is also used to
        # determine which columns to add to the general statistics table
        # https://whatshap.readthedocs.io/en/latest/guide.html#the-tsv-statistics-format
        stats_headers = {
            "variants": {
                "title": "Input Variants",
                "description": """Number of biallelic variants in the input VCF, but
                                        excluding any non-SNV variants if --only-snvs was
                                        used""",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "heterozygous_variants": {
                "title": "Heterozygous Variants",
                "description": """The number of biallelic, heterozygous variants in
                                        the input VCF. This is a subset of Input Variants.""",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "heterozygous_snvs": {
                "title": "Heterozygous SNVs",
                "description": """The number of biallelic, heterozygous SNVs in the
                                        input VCF. This is a subset of Heterozygous
                                        Variants.""",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "unphased": {
                "title": "Unphased Variants",
                "description": """The number of biallelic, heterozygous variants that
                                        are not marked as phased in the input VCF. This
                                        is a subset of heterozygous_variants.""",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "phased": {
                "title": "Phased Variants",
                "description": """The number of biallelic, heterozygous variants that
                                        are marked as phased in the input VCF. This is
                                        a subset of heterozygous_variants. Also, phased +
                                        unphased + singletons = heterozygous_variants.""",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "phased_snvs": {
                "title": "Phased SNVs",
                "description": """The number of biallelic, heterozygous SNVs that are
                                        marked as phased in the input VCF. This is a subset
                                        of phased.""",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "blocks": {
                "title": "Blocks",
                "description": "The total number of phase sets/blocks.",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "singletons": {
                "title": "Singletons",
                "description": "The number of blocks that contain exactly one variant.",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "bp_per_block_sum": {
                "title": "Total Phased bp",
                "description": """The sum of the lengths of all non-singleton
                                        blocks, where the length of a block is the
                                        number of base pairs it covers minus 1. That is, a
                                        block with two variants at positions 2 and 5 has
                                        length 3.""",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "variant_per_block_avg": {
                "title": "Avg Variants per Block",
                "description": """Description of the distribution of non-singleton
                                        block sizes, where the size of a block is the number
                                        of variants it contains. Number of biallelic
                                        variants in the input VCF, but excluding any non-SNV
                                        variants if --only-snvs was used.""",
                "hidden": False,
            },
            "frac_het_phased": {
                "title": "% Phased Variants",
                "description": """Fraction of heterozygous variants
                                          that could be phased.
                                        """,
                "format": "{:,.0f}",
                "suffix": "%",
                "min": 0,
                "max": 100,
                "modify": lambda x: x * 100,
                "hidden": True,
            },
            "bp_per_block_avg": {
                "rid": "bp_per_block_avg",
                "title": "Avg bp per Block",
                "description": """Description of the distribution of non-singleton
                                        block lengths, where the length of a block is the
                                        number of base pairs it covers minus 1. That is, a
                                        block with two variants at positions 2 and 5 has
                                        length 3.""",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "block_n50": {
                "title": "NG50",
                "description": """The NG50 value of the distribution of the block
                                        lengths. Interleaved blocks are cut in order to
                                        avoid artificially inflating this value.""",
                "format": "{:,.0f}",
                "hidden": True,
            },
        }

        general = dict()
        for sample, sample_stats in self.whatshap_stats.items():
            # Get the summary field
            summary_field = self.get_summary_field(sample_stats)

            # The summary data we are adding to the general statistics
            summary_data = sample_stats[summary_field]

            stats = {field: summary_data[field] for field in stats_headers}

            general[sample] = stats

        self.add_section(
            name="WhatsHap statistics",
            anchor="whatshap-table",
            plot=table.plot(
                general,
                stats_headers,
                {
                    "namespace": "WhatsHap",
                    "id": "whatshap-stats-table",
                    "title": "WhatsHap Statistics",
                },
            ),
        )
