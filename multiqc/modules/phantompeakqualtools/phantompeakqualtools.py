""" MultiQC module to parse output from phantompeakqualtools """


import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="phantompeakqualtools",
            anchor="phantompeakqualtools",
            href="https://www.encodeproject.org/software/phantompeakqualtools",
            info="computes informative enrichment and quality measures for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data.",
            doi=["10.1101/gr.136184.111", "10.1038/nbt.1508"],
        )

        # Parse logs
        self.phantompeakqualtools_data = dict()
        for f in self.find_log_files("phantompeakqualtools/out", filehandles=False):
            self.parse_phantompeakqualtools(f)

        # Filter to strip out ignored sample names
        self.phantompeakqualtools_data = self.ignore_samples(self.phantompeakqualtools_data)

        # Warning when no files are found
        if len(self.phantompeakqualtools_data) == 0:
            raise ModuleNoSamplesFound

        # Log
        log.info(f"Found {len(self.phantompeakqualtools_data)} logs")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None, f["s_name"])

        # Write parsed data to a file
        self.write_data_file(self.phantompeakqualtools_data, "multiqc_phantompeakqualtools")

        # Report section
        self.phantompeakqualtools_general_stats()

    # Parse spp.out file from phantompeakqualtools
    def parse_phantompeakqualtools(self, f):
        parsed_data = {}
        lines = f["f"].splitlines()
        for line in lines:
            s = line.split("\t")
            parsed_data["Estimated_Fragment_Length_bp"] = int(s[2].split(",")[0])
            parsed_data["NSC"] = float(s[8])
            parsed_data["RSC"] = float(s[9])
        if len(parsed_data) > 0:
            if f["s_name"] in self.phantompeakqualtools_data:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.add_data_source(f)
            self.phantompeakqualtools_data[f["s_name"]] = parsed_data

    # Report fragment length, NSC and RSC in general stat table
    def phantompeakqualtools_general_stats(self):
        """Add columns to General Statistics table"""
        headers = {
            "Estimated_Fragment_Length_bp": {
                "title": "Frag Length",
                "description": "Estimated fragment length (bp)",
                "min": 0,
                "format": "{:,.0f}",
            },
            "NSC": {
                "title": "NSC",
                "description": "Normalized strand cross-correlation",
                "max": 10,
                "min": 0,
                "format": "{:,.2f}",
                "scale": "RdYlGn-rev",
            },
            "RSC": {
                "title": "RSC",
                "description": "Relative strand cross-correlation",
                "max": 10,
                "min": 0,
                "format": "{:,.2f}",
                "scale": "RdYlBu-rev",
            },
        }
        self.general_stats_addcols(self.phantompeakqualtools_data, headers)
