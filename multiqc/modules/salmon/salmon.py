""" MultiQC module to parse output from Salmon """

import json
import logging
import os

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph
from multiqc.utils import config

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Salmon",
            anchor="salmon",
            href="https://combine-lab.github.io/salmon/",
            info="is a tool for quantifying the expression of transcripts using RNA-seq data.",
            doi="10.1038/nmeth.4197",
        )

        # Parse meta information. JSON win!
        self.salmon_meta = dict()
        for f in self.find_log_files("salmon/meta"):
            # Get the s_name from the parent directory
            if os.path.basename(f["root"]) in ["aux_info", "aux"]:
                s_name = os.path.basename(os.path.dirname(f["root"]))
                s_name = self.clean_s_name(s_name, f)
                self.salmon_meta[s_name] = json.loads(f["f"])
                self.add_software_version(self.salmon_meta[s_name]["salmon_version"], s_name)

        # Parse Fragment Length Distribution logs
        self.salmon_fld = dict()
        for f in self.find_log_files("salmon/fld"):
            # Get the s_name from the parent directory
            if os.path.basename(f["root"]) == "libParams":
                s_name = os.path.basename(os.path.dirname(f["root"]))
                s_name = self.clean_s_name(s_name, f)
                parsed = dict()
                for i, v in enumerate(f["f"].split()):
                    parsed[i] = float(v)
                if len(parsed) > 0:
                    if s_name in self.salmon_fld:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name)
                    self.salmon_fld[s_name] = parsed

        # Parse Library Format Counts information. JSON file expected
        self.salmon_lfc = dict()
        for f in self.find_log_files("salmon/lfc"):
            s_name = os.path.basename(f["root"])  # lfc file located at root folder
            s_name = self.clean_s_name(s_name, f)
            self.salmon_lfc[s_name] = json.loads(f["f"])

        # Filter to strip out ignored sample names
        self.salmon_meta = self.ignore_samples(self.salmon_meta)
        self.salmon_fld = self.ignore_samples(self.salmon_fld)
        self.salmon_lfc = self.ignore_samples(self.salmon_lfc)

        if len(self.salmon_meta) == 0 and len(self.salmon_fld) == 0 and len(self.salmon_lfc) == 0:
            raise ModuleNoSamplesFound

        if len(self.salmon_meta) > 0:
            log.info(f"Found {len(self.salmon_meta)} meta reports")
            self.write_data_file(self.salmon_meta, "multiqc_salmon")
        if len(self.salmon_fld) > 0:
            log.info(f"Found {len(self.salmon_fld)} fragment length distributions")
        if len(self.salmon_lfc) > 0:
            log.info(f"Found {len(self.salmon_lfc)} library format counts reports")

        if self.salmon_meta:
            # Add alignment rate to the general stats table
            # Convert library types to string:
            for d in self.salmon_meta.values():
                if "library_types" in d:
                    d["library_types"] = ", ".join(d["library_types"])

            headers = {
                "percent_mapped": {
                    "title": "% Aligned",
                    "description": "% Mapped reads",
                    "max": 100,
                    "min": 0,
                    "suffix": "%",
                    "scale": "YlGn",
                },
                "num_mapped": {
                    "title": "M Aligned",
                    "description": "Mapped reads (millions)",
                    "min": 0,
                    "scale": "PuRd",
                    "modify": lambda x: float(x) * config.read_count_multiplier,
                    "suffix": config.read_count_prefix,
                    "shared_key": "read_count",
                },
                "library_types": {
                    "title": "Library types",
                    "description": "Library types",
                    "scale": False,
                    # Hide if all samples have the same value
                    "hidden": len(set(d.get("library_types") for d in self.salmon_meta.values())) == 1,
                },
            }
            self.general_stats_addcols(self.salmon_meta, headers)

        if self.salmon_lfc:
            # Compatible fragments ratios data
            lfc_headers = {
                "compatible_fragment_ratio": {
                    "title": "CFR",
                    "description": "Compatible fragment ratio",
                    "min": 0.0,
                    "max": 100.0,
                    "modify": lambda x: x * 100.0,
                    "suffix": "%",
                    "scale": "YlGn",
                },
                "strand_mapping_bias": {
                    "title": "M Bias",
                    "description": "Strand mapping bias",
                    "scale": "BuGn",
                    "max": 1.0,
                },
            }
            # add strand mapping bias data
            self.general_stats_addcols(self.salmon_lfc, lfc_headers)

        if self.salmon_fld:
            # Fragment length distribution plot
            pconfig = {
                "smooth_points": 500,
                "id": "salmon_plot",
                "title": "Salmon: Fragment Length Distribution",
                "ylab": "Fraction",
                "xlab": "Fragment Length (bp)",
                "ymin": 0,
                "xmin": 0,
                "tt_label": "<b>{point.x:,.0f} bp</b>: {point.y:,.0f}",
            }
            self.add_section(plot=linegraph.plot(self.salmon_fld, pconfig))
