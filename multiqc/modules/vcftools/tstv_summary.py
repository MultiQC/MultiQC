"""MultiQC module to parse TsTv by summary output from vcftools TsTv-summary"""

import logging

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class TsTvSummaryMixin:
    def parse_tstv_summary(self):
        """Create the HTML for the TsTv summary plot."""

        self.vcftools_tstv_summary = dict()
        for f in self.find_log_files("vcftools/tstv_summary", filehandles=True):
            d = {}
            for line in f["f"].readlines()[1:]:  # don't add the header line (first row)
                key = line.split()[0]  # taking the first column (MODEL) as key
                val = int(line.split()[1])  # taking the second column (COUNT) as value
                d[key] = val
            self.vcftools_tstv_summary[f["s_name"]] = d
            self.add_data_source(f, "Summary")

        # Filter out ignored sample names
        self.vcftools_tstv_summary = self.ignore_samples(self.vcftools_tstv_summary)
        if len(self.vcftools_tstv_summary) == 0:
            return 0

        # Write data to file
        self.write_data_file(self.vcftools_tstv_summary, "vcftools_tstv_summary")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Specifying the categories of the bargraph
        keys = ["AC", "AG", "AT", "CG", "CT", "GT", "Ts", "Tv"]

        pconfig = {
            "id": "vcftools_tstv_summary",
            "title": "VCFTools: TsTv Summary",
            "ylab": "Counts",
        }

        self.add_section(
            name="TsTv Summary",
            anchor="vcftools-tstv-summary",
            description="Plot of `TSTV-SUMMARY` - count of different types of transition and transversion SNPs.",
            plot=bargraph.plot(self.vcftools_tstv_summary, keys, pconfig),
        )

        return len(self.vcftools_tstv_summary)
