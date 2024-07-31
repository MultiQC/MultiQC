import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="OptiType",
            anchor="optitype",
            href="https://github.com/FRED-2/OptiType",
            info="Precision HLA typing from next-generation sequencing data.",
            extra="""
            Novel HLA genotyping algorithm based on integer linear programming, capable of producing accurate 4-digit 
            HLA genotyping predictions from NGS data by simultaneously selecting all major and minor HLA Class I alleles.
            """,
            doi="10.1093/bioinformatics/btu548",
        )

        # Find and load any optitype reports
        self.optitype_data = dict()

        for f in self.find_log_files("optitype"):
            rows = f["f"].splitlines()
            # First col is empty / or always zero
            headers = rows[0].split("\t")[1:]
            cols = rows[1].split("\t")[1:]

            if f["s_name"] in self.optitype_data:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.optitype_data[f["s_name"]] = dict()
            for header, col in zip(headers, cols):
                self.optitype_data[f["s_name"]][header] = col
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.optitype_data = self.ignore_samples(self.optitype_data)

        if len(self.optitype_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.optitype_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.optitype_data, "multiqc_optitype")

        # Create OptiType overview table
        self.add_summary_metrics()

        # Create bar graph of subtype sample counts
        self.summary_barplot()

    def add_summary_metrics(self):
        """Take the parsed entries from OptiType and add them to the main plot"""
        # A1, A2, B1, B2, C1, C2, Reads, Objective
        headers = {
            "A1": {
                "title": "HLA-A1",
                "description": "First HLA-A allele",
                "scale": False,
            },
            "A2": {
                "title": "HLA-A2",
                "description": "Second HLA-A allele",
                "scale": False,
                "hidden": True,
            },
            "B1": {
                "title": "HLA-B1",
                "description": "First HLA-B allele",
                "scale": False,
            },
            "B2": {
                "title": "HLA-B2",
                "description": "Second HLA-B allele",
                "scale": False,
                "hidden": True,
            },
            "C1": {
                "title": "HLA-C1",
                "description": "First HLA-C allele",
                "scale": False,
            },
            "C2": {
                "title": "HLA-C2",
                "description": "Second HLA-C allele",
                "scale": False,
                "hidden": True,
            },
            "Reads": {
                "title": "Reads",
                "description": "Number of reads covering the HLA",
                "scale": "YlGnBu",
                "hidden": True,
                "format": "{:,.0f}",
            },
            "Objective": {
                "title": "Objective Score",
                "description": "Score of the objective function for the prediction.",
                "scale": "BuGn",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.optitype_data, headers)

    def summary_barplot(self):
        """Make a bar plot showing the number of samples for each allele"""

        alleles = ["A1", "A2", "B1", "B2", "C1", "C2"]
        data = {al: {} for al in alleles}
        for s_data in self.optitype_data.values():
            for al in alleles:
                if s_data[al] not in data[al]:
                    data[al][s_data[al]] = 0
                data[al][s_data[al]] += 1

        pconfig = {
            "id": "optitype_summary_plot",
            "title": "OptiType: Summary of alleles",
            "cpswitch": False,
            "ylab": "# Samples",
            "y_decimals": False,
            "use_legend": False,
        }

        self.add_section(
            name="Summary of alleles",
            description="Number of samples sharing the same allele.",
            anchor="optitype_summary",
            plot=bargraph.plot(data, pconfig=pconfig),
        )
