#!/usr/bin/env python

""" MultiQC module to parse output from bustools inspect """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import json

from multiqc.plots import bargraph, table
from multiqc.utils import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Bustools",
            anchor="bustools",
            href="https://bustools.github.io/",
            info="works with BUS files - a file format for single-cell RNA-seq data "
            "designed to facilitate the development of modular workflows for data "
            "processing.",
            doi="10.1093/bioinformatics/btz279",
        )

        self.prepare_data()

        log.info("Found {} logs".format(len(self.bustools_data)))
        self.write_data_file(self.bustools_data, "multiqc_macs")

        self.bustools_general_stats()
        self.bustools_section()

    def prepare_data(self):
        # Parse logs
        self.bustools_data = dict()
        for f in self.find_log_files("bustools", filehandles=True):
            content = json.load(f["f"])
            s_name = self.clean_s_name(os.path.basename(f["root"]), f, root=os.path.dirname(f["root"]))
            self.bustools_data[s_name] = content
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.bustools_data = self.ignore_samples(self.bustools_data)

        if len(self.bustools_data) == 0:
            raise UserWarning

        # now fill out the table(s) headers
        self.headers = OrderedDict()
        self.headers["numRecords"] = {
            "title": "{} Bus Records".format(config.read_count_prefix),
            "description": "The number of Bus Records ({})".format(config.read_count_desc),
            "scale": "BuGn",
            "min": 0,
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "format": "{:,.2f}",
        }
        self.headers["numReads"] = {
            "title": "{} Reads".format(config.read_count_prefix),
            "description": "Total number of reads ({})".format(config.read_count_desc),
            "min": 0,
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "format": "{:,.2f}",
        }
        self.headers["numBarcodes"] = {
            "title": "Barcodes",
            "description": "Number of distinct barcodes",
            "min": 0,
            "scale": "YlGn",
            "format": "{:,.0f}",
            "shared_key": "barcodes",
        }
        self.headers["meanReadsPerBarcode"] = {
            "title": "Mean reads per barcode",
            "scale": "BuGn",
            "min": 0,
            "format": "{:,.2f}",
        }
        self.headers["numUMIs"] = {
            "title": "{} Distinct UMIs".format(config.read_count_prefix),
            "description": "Number of distinct Unique Molecular Identifiers (UMIs) ({})".format(config.read_count_desc),
            "scale": "Purples",
            "min": 0,
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "format": "{:,.2f}",
        }
        self.headers["numBarcodeUMIs"] = {
            "title": "{} Distinct barcode-UMI".format(config.read_count_prefix),
            "description": "Number of distinct barcode and Unique Molecular Identifiers (UMIs) pairs ({})".format(
                config.read_count_desc
            ),
            "scale": "Greens",
            "min": 0,
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "format": "{:,.2f}",
        }
        self.headers["meanUMIsPerBarcode"] = {
            "title": "Mean UMIs per barcode",
            "scale": "PuBnGn",
            "min": 0,
            "format": "{:,.2f}",
        }
        self.headers["gtRecords"] = {
            "title": "{} 2xdepth records".format(config.read_count_prefix),
            "description": "Estimated number of new records at 2x sequencing depth ({})".format(config.read_count_desc),
            "min": 0,
            "scale": "Oranges",
            "format": "{:,.2f}",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.headers["numTargets"] = {
            "title": "{} Distinct targets".format(config.read_count_prefix),
            "description": "Number of distinct targets detected ({})".format(config.read_count_desc),
            "min": 0,
            "scale": "BuGn",
            "hidden": True,
            "format": "{:,.2f}",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.headers["meanTargetsPerSet"] = {
            "title": "Mean targets",
            "description": "Mean number of targets per set",
            "min": 0,
            "scale": "Greens",
            "hidden": True,
            "format": "{:,.2f}",
        }
        self.headers["numSingleton"] = {
            "title": "{} Singleton targets".format(config.read_count_prefix),
            "description": "Number of reads with singleton target ({})".format(config.read_count_desc),
            "min": 0,
            "scale": "Blues",
            "hidden": True,
            "format": "{:,.2f}",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.headers["gtTargets"] = {
            "title": "{} 2xdepth targets".format(config.read_count_prefix),
            "description": "Estimated number of new targets at 2x sequencing depth ({})".format(config.read_count_desc),
            "min": 0,
            "scale": "BuGn",
            "hidden": True,
            "format": "{:,.2f}",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.headers["numBarcodesOnWhitelist"] = {
            "title": "Whitelisted barcodes",
            "description": "Number of barcodes in agreement with whitelist",
            "min": 0,
            "scale": "Greens",
            "hidden": True,
            "format": "{:,.0f}",
            "shared_key": "barcodes",
        }
        self.headers["percentageBarcodesOnWhitelist"] = {
            "title": "% Whitelisted barcodes",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
        }
        self.headers["numReadsOnWhitelist"] = {
            "title": "{} Whitelisted reads".format(config.read_count_prefix),
            "description": "Number of reads with barcode in agreement with whitelist ({})".format(
                config.read_count_desc
            ),
            "scale": "PuBu",
            "min": 0,
            "hidden": True,
            "format": "{:,.2f}",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.headers["percentageReadsOnWhitelist"] = {
            "title": "% Whitelisted reads",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
        }

    def bustools_general_stats(self):
        # we take a smaller subset
        general_stats_headers = [
            "numRecords",
            "numBarcodes",
            "meanReadsPerBarcode",
            "numUMIs",
            "gtRecords",
            "percentageReadsOnWhitelist",
        ]
        headers = OrderedDict()
        for header, v in self.headers.items():
            if header in general_stats_headers:
                headers[header] = {key: value for key, value in v.items()}  # deep copy
                if header not in ["numRecords", "numUMIs", "percentageReadsOnWhitelist"]:
                    headers[header]["hidden"] = True
        self.general_stats_addcols(self.bustools_data, headers)

    def bustools_section(self):
        """Add bargraphs showing the mean UMIs per barcode and percentages in whitelist"""
        # add the summary table
        tconfig = {"namespace": "Bustools", "id": "bustools_summary", "table_title": "Bustools Summary Table"}
        self.add_section(
            name="Summary table",
            anchor="bustools-inspect",
            description="This is a table of the complete output of bustools inspect. Note that some columns are hidden by default (click <em>Configure Columns</em> to show).",
            plot=table.plot(self.bustools_data, self.headers, tconfig),
        )

        # also make some nice barplots
        # barplot for mean umis per sample
        mean_umis = {
            sample: {"UMIs per barcode": values["meanUMIsPerBarcode"]} for sample, values in self.bustools_data.items()
        }

        self.add_section(
            name="Mean number of UMIs per barcode",
            anchor="bustools-umis",
            description="Average number of UMIs (unique molecular identifiers) per barcode",
            helptext="Each unique barcode represents a cell and each Unique Molecular Identifier (UMI) represents "
            "a unique transcript molecule. By counting the mean number of UMIs per barcode, you "
            "effectively calculate the average number of unique transcripts per cell.",
            plot=bargraph.plot(
                mean_umis,
                pconfig={
                    "id": "bus_umis",
                    "title": "Bustools: Mean number of UMIs per barcode per sample",
                    "cpswitch": False,
                    "tt_percentages": False,
                    "ylab": "Mean UMIs per barcode",
                },
            ),
        )

        # barplot for the percentage of reads and barcodes on the whitelist
        percentage_whitelist = {
            sample: {
                "Reads on whitelist": values["percentageReadsOnWhitelist"],
                "Barcodes on whitelist": values["percentageBarcodesOnWhitelist"],
            }
            for sample, values in self.bustools_data.items()
        }
        self.add_section(
            name="Percentage in whitelist",
            anchor="bustools-reads",
            description="The whitelist is a list of unique barcodes used in your protocol, either provided or inferred from the data.",
            helptext="Each unique barcode from the whitelist represents a cell. The percentage of "
            "reads with barcode / barcodes in the whitelist is a measure of percentage of reads that could "
            "be asigned to a cell.",
            plot=bargraph.plot(
                percentage_whitelist,
                pconfig={
                    "id": "bus_reads",
                    "title": "Bustools: Barcodes / reads with barcodes in the whitelist",
                    "ymax": 100,
                    "ymix": 0,
                    "cpswitch": False,
                    "tt_percentages": False,
                    "ylab": "Percentage of barcodes / reads with barcodes in the whitelist",
                    "stacking": None,
                    "ylab_format": "{value}%",
                },
            ),
        )
