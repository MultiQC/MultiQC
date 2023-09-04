""" MultiQC module to parse output from WhatsHap """

import logging
from collections import OrderedDict, defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    WhatsHap module class, parses WhatsHap output.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="WhatsHap",
            anchor="whatshap",
            href="https://whatshap.readthedocs.io/",
            info="""is a program for phasing genomic variants using DNA sequencing
                    reads, also called read-based phasing or haplotype assembly. It
                    is especially suitable for long reads, but also works well with
                    short reads.
                    """,
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

        # Raise UserWarning if we didn't find any data
        if not self.whatshap_stats:
            raise UserWarning

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

        def process_data(data):
            """
            Process the data to convert it to the proper format

            In practice, this only involves the following changes:
                - Convert all numeric values to integers
                - Replace 'nan' values with zero
            """
            for key, value in data.items():
                # Try to make the value an int
                try:
                    data[key] = int(data[key])
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
            spline = line.strip().split("\t")
            data = {field: value for field, value in zip(header, spline)}

            # Remove the sample and chromsome from the data
            sample = str(data.pop("sample"))
            chromosome = data.pop("chromosome")

            # Process the remaining data fields
            process_data(data)

            # Calculate the fraction of heterozygous variants that were phased
            try:
                frac_het_phased = data["phased"] / data["heterozygous_variants"]
            except ZeroDivisionError:
                frac_het_phased = 0
            data["frac_het_phased"] = frac_het_phased

            # Insert the current line under chromosome
            results[chromosome] = data

        # If we were parsing an empty file, we return dummy data (basically all
        # zero) so the sample statistics are consistent, which eases parsing
        # later on. WhatsHap uses the fictional 'ALL' chromsome to store the
        # summary for all results, so we re-use that to store all zero's
        if not results:
            results = dict()
            results["ALL"] = defaultdict(int)

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
            summary_field = "ALL"
        return summary_field

    def whatshap_add_general_stats(self):
        """Add WhatsHap stats to the general statistics table"""

        # Store the configuration for each column, this is also used to
        # determine which columns to add to the general statistics table
        # https://whatshap.readthedocs.io/en/latest/guide.html#the-tsv-statistics-format
        general_stats_headers = OrderedDict(
            [
                (
                    "frac_het_phased",
                    {
                        "id": "perc_het_phased",
                        "title": "% Phased Variants",
                        "description": """Fraction of heterozygous variants
                                          that could be phased.
                                        """,
                        "format": "{:,.0f}",
                        "suffix": "%",
                        "min": 0,
                        "max": 100,
                        "modify": lambda x: x * 100,
                        "hidden": False,
                    },
                ),
                (
                    "bp_per_block_avg",
                    {
                        "id": "bp_per_block_avg",
                        "title": "Avg bp per Block",
                        "description": """Description of the distribution of non-singleton
                                        block lengths, where the length of a block is the
                                        number of base pairs it covers minus 1. That is, a
                                        block with two variants at positions 2 and 5 has
                                        length 3.""",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "block_n50",
                    {
                        "id": "block_n50",
                        "title": "NG50",
                        "description": """The NG50 value of the distribution of the block
                                        lengths. Interleaved blocks are cut in order to
                                        avoid artificially inflating this value.""",
                        "format": "{:,.0f}",
                        "hidden": True,
                    },
                ),
            ]
        )

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
            "anchor": "multiqc_whatshap_phased_bp",
            "ylab": "Base Pairs",
            "cpswitch": False,
            "tt_percentages": False,
        }

        keys = OrderedDict()
        keys["Phased Base Pairs"] = {"name": "Phased Base Pairs"}

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
        stats_headers = OrderedDict(
            [
                (
                    "variants",
                    {
                        "id": "variants",
                        "title": "Input Variants",
                        "description": """Number of biallelic variants in the input VCF, but
                                        excluding any non-SNV variants if --only-snvs was
                                        used""",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "heterozygous_variants",
                    {
                        "id": "heterozygous_variants",
                        "title": "Heterozygous Variants",
                        "description": """The number of biallelic, heterozygous variants in
                                        the input VCF. This is a subset of Input Variants.""",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "heterozygous_snvs",
                    {
                        "id": "heterozygous_snvs",
                        "title": "Heterozygous SNVs",
                        "description": """The number of biallelic, heterozygous SNVs in the
                                        input VCF. This is a subset of Heterozygous
                                        Variants.""",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "unphased",
                    {
                        "id": "unphased",
                        "title": "Unphased Variants",
                        "description": """The number of biallelic, heterozygous variants that
                                        are not marked as phased in the input VCF. This
                                        is a subset of heterozygous_variants.""",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "phased",
                    {
                        "id": "phased",
                        "title": "Phased Variants",
                        "description": """The number of biallelic, heterozygous variants that
                                        are marked as phased in the input VCF. This is
                                        a subset of heterozygous_variants. Also, phased +
                                        unphased + singletons = heterozygous_variants.""",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "phased_snvs",
                    {
                        "id": "phased_snvs",
                        "title": "Phased SNVs",
                        "description": """The number of biallelic, heterozygous SNVs that are
                                        marked as phased in the input VCF. This is a subset
                                        of phased.""",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "blocks",
                    {
                        "id": "blocks",
                        "title": "Blocks",
                        "description": "The total number of phase sets/blocks.",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "singletons",
                    {
                        "id": "singletons",
                        "title": "Singletons",
                        "description": "The number of blocks that contain exactly one variant.",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "bp_per_block_sum",
                    {
                        "id": "bp_per_block_sum",
                        "title": "Total Phased bp",
                        "description": """The sum of the lengths of all non-singleton
                                        blocks, where the length of a block is the
                                        number of base pairs it covers minus 1. That is, a
                                        block with two variants at positions 2 and 5 has
                                        length 3.""",
                        "format": "{:,.0f}",
                        "hidden": False,
                    },
                ),
                (
                    "variant_per_block_avg",
                    {
                        "id": "variant_per_block_avg",
                        "title": "Avg Variants per Block",
                        "description": """Description of the distribution of non-singleton
                                        block sizes, where the size of a block is the number
                                        of variants it contains. Number of biallelic
                                        variants in the input VCF, but excluding any non-SNV
                                        variants if --only-snvs was used.""",
                        "hidden": False,
                    },
                ),
                (
                    "frac_het_phased",
                    {
                        "id": "perc_het_phased",
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
                ),
                (
                    "bp_per_block_avg",
                    {
                        "id": "bp_per_block_avg",
                        "title": "Avg bp per Block",
                        "description": """Description of the distribution of non-singleton
                                        block lengths, where the length of a block is the
                                        number of base pairs it covers minus 1. That is, a
                                        block with two variants at positions 2 and 5 has
                                        length 3.""",
                        "format": "{:,.0f}",
                        "hidden": True,
                    },
                ),
                (
                    "block_n50",
                    {
                        "id": "block_n50",
                        "title": "NG50",
                        "description": """The NG50 value of the distribution of the block
                                        lengths. Interleaved blocks are cut in order to
                                        avoid artificially inflating this value.""",
                        "format": "{:,.0f}",
                        "hidden": True,
                    },
                ),
            ]
        )

        general = dict()
        for sample, sample_stats in self.whatshap_stats.items():
            # Get the summary field
            summary_field = self.get_summary_field(sample_stats)

            # The summary data we are adding to the general statistics
            summary_data = sample_stats[summary_field]

            stats = {field: summary_data[field] for field in stats_headers}

            general[sample] = stats

        self.add_section(name="WhatsHap statistics", anchor="whatshap-table", plot=table.plot(general, stats_headers))
