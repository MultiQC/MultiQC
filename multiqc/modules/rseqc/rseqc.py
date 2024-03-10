""" MultiQC module to parse output from RSeQC """

import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """RSeQC is a collection of scripts. This MultiQC module
    supports some but not all. The code for each script is split
    into its own file and adds a section to the module output if
    logs are found."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="RSeQC",
            anchor="rseqc",
            href="http://rseqc.sourceforge.net/",
            info="package provides a number of useful modules that can"
            " comprehensively evaluate high throughput RNA-seq data.",
            doi="10.1093/bioinformatics/bts356",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # Get the list of submodules (can be customised)
        rseqc_sections = getattr(config, "rseqc_sections", [])
        if len(rseqc_sections) == 0:
            rseqc_sections = [
                "read_distribution",
                "gene_body_coverage",
                "inner_distance",
                "read_gc",
                "read_duplication",
                "junction_annotation",
                "junction_saturation",
                "infer_experiment",
                "bam_stat",
                "tin",
            ]

        # Call submodule functions
        for sm in rseqc_sections:
            try:
                # Import the submodule and call parse_reports()
                #   Function returns number of parsed logs
                module = __import__(f"multiqc.modules.rseqc.{sm}", fromlist=[""])
                n[sm] = getattr(module, "parse_reports")(self)
                if n[sm] > 0:
                    log.info(f"Found {n[sm]} {sm} reports")
            except (ImportError, AttributeError):
                log.error(f"Could not find RSeQC Section '{sm}'")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        if max((len(vals) for vals in self.general_stats_data.values()), default=0) > 0:
            self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
