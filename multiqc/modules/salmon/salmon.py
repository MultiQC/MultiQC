#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function
from collections import OrderedDict
import json
import logging
import os

from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Salmon",
            anchor="salmon",
            href="http://combine-lab.github.io/salmon/",
            info="is a tool for quantifying the expression of transcripts using RNA-seq data.",
        )

        # Parse meta information. JSON win!
        self.salmon_meta = dict()
        for f in self.find_log_files("salmon/meta"):
            # Get the s_name from the parent directory
            s_name = os.path.basename(os.path.dirname(f["root"]))
            s_name = self.clean_s_name(s_name, f)
            self.salmon_meta[s_name] = json.loads(f["f"])

        # Parse Fragment Length Distribution logs
        self.salmon_fld = dict()
        for f in self.find_log_files("salmon/fld"):
            # Get the s_name from the parent directory
            if os.path.basename(f["root"]) == "libParams":
                s_name = os.path.basename(os.path.dirname(f["root"]))
                s_name = self.clean_s_name(s_name, f)
                parsed = OrderedDict()
                for i, v in enumerate(f["f"].split()):
                    parsed[i] = float(v)
                if len(parsed) > 0:
                    if s_name in self.salmon_fld:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.salmon_fld[s_name] = parsed

        # Parse Library Format Counts information. JSON file expected
        self.salmon_lfc = dict()
        for f in self.find_log_files("salmon"):
            s_name = os.path.basename(f["root"])  # lfc file located at root folder
            s_name = self.clean_s_name(s_name, f)
            self.salmon_lfc[s_name] = json.loads(f["f"])

        # Filter to strip out ignored sample names
        self.salmon_meta = self.ignore_samples(self.salmon_meta)
        self.salmon_fld = self.ignore_samples(self.salmon_fld)
        self.salmon_lfc = self.ignore_samples(self.salmon_lfc)

        if len(self.salmon_meta) == 0 and len(self.salmon_fld) == 0:
            raise UserWarning

        if len(self.salmon_meta) > 0:
            log.info("Found {} meta reports".format(len(self.salmon_meta)))
            self.write_data_file(self.salmon_meta, "multiqc_salmon")
        if len(self.salmon_fld) > 0:
            log.info("Found {} fragment length distributions".format(len(self.salmon_fld)))

        # Add alignment rate to the general stats table
        headers = OrderedDict()
        headers["percent_mapped"] = {
            "title": "% Aligned",
            "description": "% Mapped reads",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlGn",
        }
        headers["num_mapped"] = {
            "title": "M Aligned",
            "description": "Mapped reads (millions)",
            "min": 0,
            "scale": "PuRd",
            "modify": lambda x: float(x) / 1000000,
            "shared_key": "read_count",
        }
        # reformat library types and check if all samples have the same value
        lib_types = None
        hide_lib = True
        for s in self.salmon_meta:
            self.salmon_meta[s]["library_types"] = ",".join(self.salmon_meta[s]["library_types"])
            if not lib_types:
                lib_types = self.salmon_meta[s]["library_types"]
            if hide_lib and self.salmon_meta[s]["library_types"] != lib_types:
                hide_lib = False
        headers["library_types"] = {
            "title": "Lib Types",
            "description": "Library types",
            "scale": False,
            "hidden": hide_lib,
        }
        self.general_stats_addcols(self.salmon_meta, headers)

        # add compatible fragments ratios data
        lfc_headers = OrderedDict()
        lfc_headers["compatible_fragment_ratio"] = {
            "title": "CFR",
            "description": "Compatible fragment ratio",
            "min": 0.0,
            "max": 1.0,
            "scale": "YlGn",
        }
        # add strand mapping bias data
        lfc_headers["strand_mapping_bias"] = {
            "title": "M Bias",
            "description": "Strand mapping bias",
            "scale": "BuGn",
            "max": 1.0,
        }
        self.general_stats_addcols(self.salmon_lfc, lfc_headers)

        if len(self.salmon_fld) > 0:
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
