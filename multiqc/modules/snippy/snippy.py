import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The following commands are implemented:

    - `snippy`
      - Variant type descriptive statistics.
      - Parses summary `prefix.txt` file that is generated.
    - `snippy-core`
      - Core genome alignment descriptive statistics.
      - Parses summary `prefix.txt` file that is generated.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Snippy",
            anchor="snippy",
            href="https://github.com/tseemann/snippy",
            info="Rapid haploid variant calling and core genome alignment.",
            # Can't find a DOI // doi=
        )

        data_by_sample: Dict[str, Dict] = {}
        core_data_by_sample: Dict[str, Dict] = {}

        self.snippy_col = [
            "Variant-COMPLEX",
            "Variant-DEL",
            "Variant-INS",
            "Variant-SNP",
            "VariantTotal",
        ]
        self.snippy_core_col = [
            "LENGTH",
            "ALIGNED",
            "UNALIGNED",
            "VARIANT",
            "HET",
            "MASKED",
            "LOWCOV",
        ]

        # Parse the txt files from snippy
        for f in self.find_log_files("snippy/snippy"):
            # Check for duplicate sample names
            if f["s_name"] in data_by_sample:
                log.debug(f"Duplicate sample name found for snippy! Overwriting: {f['s_name']}")
            # Add the file data under the key filename
            data = self.parse_snippy_txt(f["f"])
            if data:
                data_by_sample[f["s_name"]] = data
                self.add_software_version(data["version"], f["s_name"])

            self.add_data_source(f, section="snippy")

        # Ignore samples specified by the user
        data_by_sample = self.ignore_samples(data_by_sample)

        # Parse the txt files from snippy-core
        for f in self.find_log_files("snippy/snippy-core"):
            # Check for duplicate sample names
            if f["s_name"] in core_data_by_sample:
                log.debug(f"Duplicate sample name found for snippy-core! Overwriting: {f['s_name']}")
            # Add the file data under the key filename
            core_data_by_sample[f["s_name"]] = self.parse_snippy_core_txt(f["f"])

            self.add_data_source(f, section="snippy-core")

        # Ignore samples specified by the user
        core_data_by_sample = self.ignore_samples(core_data_by_sample)

        # Raise warning if no logs were found.
        if len(data_by_sample) == 0 and len(core_data_by_sample) == 0:
            raise ModuleNoSamplesFound

        # Run analysis if txt files found
        if len(data_by_sample) > 0:
            log.info(f"Found {len(data_by_sample)} reports")
            self.snippy_stats_table(data_by_sample)
            self.snippy_report_section(data_by_sample)

        if len(core_data_by_sample) > 0:
            log.info(f"Found {len(core_data_by_sample)} snippy-core reports")
            self.snippy_core_stats_table(core_data_by_sample)
            self.snippy_core_report_section(core_data_by_sample)

    def parse_snippy_txt(self, file):
        """
        Parse the txt file for snippy output.
        """
        data = {}
        for line in file.splitlines():
            split_line = line.strip().split("\t")
            if split_line[0] in self.snippy_col:
                data[split_line[0]] = int(split_line[1])

            if split_line[0] == "Software":
                data["version"] = split_line[1].split(" ")[1]

        if len(data) == 0:
            return False
        for col in self.snippy_col:
            if col not in data:
                data[col] = 0
        return data

    def parse_snippy_core_txt(self, file):
        """
        Parse the txt file for snippy-core output.
        """
        data = {}
        for line in file.splitlines()[1:]:
            split_line = line.split("\t")
            data[split_line[0]] = split_line[1:]
        return data

    def snippy_stats_table(self, data_by_sample):
        """
        Add the snippy data to the general stats.
        """
        self.general_stats_addcols(
            data_by_sample,
            {
                "VariantTotal": {
                    "title": "# Variants",
                    "description": "Total variants detected.",
                    "scale": "BuPu",
                    "format": "{:.0f}",
                    "shared_key": "variant_count",
                }
            },
        )
        self.write_data_file(data_by_sample, "multiqc_snippy")

    def snippy_core_stats_table(self, data_by_sample):
        """
        Add the snippy-core data to the general stats.
        """
        # Parse the statistics
        data: Dict = {}
        for file in data_by_sample:
            for sample in data_by_sample[file]:
                data[sample] = {}
                for i in range(0, len(self.snippy_core_col)):
                    data[sample][self.snippy_core_col[i]] = int(data_by_sample[file][sample][i])
                data[sample]["Percent_Aligned"] = (data[sample]["ALIGNED"] / data[sample]["LENGTH"]) * 100
                data[sample]["Percent_Het"] = (data[sample]["HET"] / data[sample]["ALIGNED"]) * 100

        self.general_stats_addcols(data, self.snippy_core_headers_config())
        self.write_data_file(data, "multiqc_snippy")

    def snippy_report_section(self, data_by_sample):
        """
        Create a report section for the snippy data.
        """
        bargraph_data: Dict[str, Dict] = {}
        for sample in data_by_sample:
            # Remove the VariantTotal stat
            bargraph_data[sample] = {}
            for stat in data_by_sample[sample]:
                if stat == "VariantTotal":
                    continue
                bargraph_data[sample][stat] = data_by_sample[sample][stat]

        # Config for the plot
        pconfig = {
            "id": "snippy_variants",
            "title": "Snippy: Variants Counts",
            "ylab": "# Variants",
        }

        self.add_section(
            name="Snippy Variants",
            anchor="snippy_variants_stats",
            plot=bargraph.plot(bargraph_data, pconfig=pconfig),
            description="This plot shows the different variant types reported by snippy.",
        )

    def snippy_core_report_section(self, data_by_sample):
        """
        Create a report section for the snippy-core data.
        """
        bargraph_data: Dict[str, Dict] = {}
        for file in data_by_sample:
            for sample in data_by_sample[file]:
                bargraph_data[sample] = {}
                for i in range(0, len(self.snippy_core_col)):
                    if self.snippy_core_col[i] == "LENGTH":
                        continue
                    bargraph_data[sample][self.snippy_core_col[i]] = int(data_by_sample[file][sample][i])
                # Subtract variant from aligned
                bargraph_data[sample]["ALIGNED"] = bargraph_data[sample]["ALIGNED"] - bargraph_data[sample]["VARIANT"]

        # Config for the plot
        pconfig = {
            "id": "snippy_core_alignment",
            "title": "Snippy: Core Alignment Statistics",
            "ylab": "# bp",
        }

        self.add_section(
            name="Snippy-Core Alignment Statistics",
            anchor="snippy_core_alignment_stats",
            plot=bargraph.plot(bargraph_data, pconfig=pconfig),
            description="Different categories of sites in a snippy-core genome alignment.",
        )

    @staticmethod
    def snippy_core_headers_config():
        """
        Prepare the headers for the snippy core stats table.
        """
        # General stats table headers
        headers = {
            "Percent_Het": {
                "title": "% Het",
                "description": "Percent of aligned sites that are heterozygous",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn-rev",
                "suffix": "%",
            },
            "HET": {
                "title": "# Hets",
                "description": "Number of heterozygous sites",
                "scale": "RdYlGn-rev",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "VARIANT": {
                "title": "# Variants",
                "description": "Number of variants",
                "scale": "Blues",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "Percent_Aligned": {
                "title": "% Aligned",
                "description": "Percent of bases aligned to the reference",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn",
                "suffix": "%",
            },
            "ALIGNED": {
                "title": "# Aligned",
                "description": "Number of aligned nucleotides",
                "scale": "RdYlGn",
                "format": "{:,.0f}",
                "hidden": True,
            },
        }
        return headers
