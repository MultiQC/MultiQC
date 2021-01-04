#!/usr/bin/env python

""" MultiQC module to parse output from ngsderive """

from __future__ import print_function
import csv
import io
import logging
import os
import re

from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, heatmap, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    ngsderive module class, parses stderr logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="ngsderive",
            anchor="ngsderive",
            href="https://github.com/stjudecloud/ngsderive",
            info="attempts to predict library information from next-generation sequencing data.",
        )

        self.strandedness = {}
        self.instrument = {}
        self.readlen = {}

        # parse ngsderive summary file
        expected_header_count_strandedness = 5
        expected_header_count_instrument = 4
        expected_header_count_readlen = 4

        for f in self.find_log_files("ngsderive/strandedness"):
            self.parse(
                self.strandedness,
                f,
                "strandedness",
                expected_header_count_strandedness,
            )

        for f in self.find_log_files("ngsderive/instrument"):
            self.parse(
                self.instrument,
                f,
                "instrument",
                expected_header_count_instrument,
            )

        for f in self.find_log_files("ngsderive/readlen"):
            self.parse(
                self.readlen,
                f,
                "readlen",
                expected_header_count_readlen,
            )

        self.strandedness = self.ignore_samples(self.strandedness)
        self.instrument = self.ignore_samples(self.instrument)
        self.readlen = self.ignore_samples(self.readlen)

        any_results_found = False
        if self.strandedness:
            self.add_strandedness_data()
            any_results_found = True

        if self.instrument:
            self.add_instrument_data()
            any_results_found = True

        if self.readlen:
            self.add_readlen_data()
            any_results_found = True

        if not any_results_found:
            raise UserWarning

    def probe_file_for_dictreader_kwargs(self, f, expected_header_count):
        """In short, this function was created to figure out which
        kwargs need to be passed to csv.DictReader. First,
        it tries to extract commonly known delimiters without
        doing an expensive scanning operation. If that
        doesn't work, we try sniffing for a dialect."""

        header = f.readline()
        f.seek(0)

        for delim in ["\t", ",", " ", "|"]:
            if delim in header and len(header.split(delim)) == expected_header_count:
                return {"delimiter": delim}

        log.warn("Could not easily detect delimiter for file. Trying csv.Sniffer()")
        dialect = csv.Sniffer().sniff(f.read(4096))  # 4096 is the buffer size
        f.seek(0)
        return {"dialect": dialect}

    def parse(self, sample_dict, found_file, subcommand, expected_header_count):
        kwargs = self.probe_file_for_dictreader_kwargs(io.StringIO(found_file["f"]), expected_header_count)
        relevant_items = []
        for row in csv.DictReader(io.StringIO(found_file["f"]), **kwargs):
            if not row.get("File"):
                continue
            sample_name = self.clean_s_name(row.get("File"), found_file["root"])
            if sample_name in sample_dict:
                log.debug("Duplicate sample name found! Overwriting: {}".format(sample_name))

            sample_dict[sample_name] = row
            self.add_data_source(f=found_file, s_name=sample_name)

    def add_strandedness_data(self):
        data = {}
        for sample, strandedness in self.strandedness.items():
            data[sample] = {
                "predicted": strandedness.get("Predicted"),
                "forward": round(float(strandedness.get("ForwardPct")) * 100.0, 2),
                "reverse": round(float(strandedness.get("ReversePct")) * 100.0, 2),
            }

        bardata = OrderedDict()
        sorted_data = sorted(data.items(), key=lambda x: x[1].get("forward"))
        for (k, v) in sorted_data:
            bardata[k] = v

        headers = OrderedDict()
        headers["predicted"] = {
            "title": "Strandedness",
            "description": "Predicted strandedness from ngsderive",
        }
        headers["forward"] = {
            "title": "% Forward Strandedness",
            "description": "Percentage of reads which were evidence for forward-stranded.",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
            "hidden": True,
        }
        headers["reverse"] = {
            "title": "% Reverse Strandedness",
            "description": "Percentage of reads which were evidence for reverse-stranded.",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
            "hidden": True,
        }
        self.general_stats_addcols(data, headers)

        # Config for the plot
        pconfig = {
            "id": "ngsderive_strandedness_plot",
            "title": "ngsderive: Strandedness",
            "ylab": "% Read Evidence",
            "ymin": 0,
            "ymax": 100,
            "tt_percentages": True,
            "ylab_format": "{value}%",
            "cpswitch": False,
        }

        self.add_section(
            name="Strandedness",
            anchor="ngsderive-strandedness",
            description="""Predicted strandedness provided by ngsderive. For more information, please see
            [the documentation](https://stjudecloud.github.io/ngsderive/subcommands/strandedness/).""",
            plot=bargraph.plot(bardata, headers, pconfig),
        )

    def add_instrument_data(self):
        general_data = {}
        for sample, instrument_data in self.instrument.items():
            general_data[sample] = {
                "instrument": " / ".join(sorted(instrument_data.get("Instrument").split(" or "))),
                "confidence": instrument_data.get("Confidence"),
                "basis": instrument_data.get("Basis"),
            }

        general_headers = OrderedDict()
        general_headers["instrument"] = {
            "title": "Predicted Instrument",
            "description": "Predicted instrument from ngsderive",
        }
        general_headers["confidence"] = {
            "title": "Instrument: Confidence",
            "description": "Level of confidence (low, medium, high) that the predicted instrument is correct.",
        }
        general_headers["basis"] = {
            "title": "Instrument: Basis",
            "description": "Basis upon which the prediction was made.",
            "hidden": True,
        }
        self.general_stats_addcols(general_data, general_headers)

        samples = []
        instruments = set()

        for s, d in general_data.items():
            samples.append(s)
            instruments.update(d.get("instrument").split(" / "))

        # move multiple instruments to the end if it exists
        instruments = sorted(instruments)
        if "multiple instruments" in instruments:
            instruments.remove("multiple instruments")
            instruments.append("multiple instruments")

        headers = OrderedDict()
        for instrument in instruments:
            headers[instrument] = {
                "title": instrument,
                "description": f"Predicted {instrument} from ngsderive",
            }

        table_data = {}
        for sample, instrument_data in self.instrument.items():
            table_data[sample] = {}
            for instrument in instrument_data.get("Instrument").split(" or "):
                table_data[sample][instrument] = instrument_data.get("Confidence")

        # Config for the plot
        config = {
            "id": "ngsderive_instruments_plot",
            "title": "ngsderive: Instruments",
        }

        self.add_section(
            name="Instrument",
            anchor="ngsderive-instrument",
            description="""Predicted instrument provided by ngsderive. For more information, please see
            [the documentation](https://stjudecloud.github.io/ngsderive/subcommands/instrument/).""",
            plot=table.plot(table_data, headers, config),
        )

    def add_readlen_data(self):
        data = {}
        for sample, readlen in self.readlen.items():
            data[sample] = {
                "evidence": readlen.get("Evidence"),
                "majoritypctdetected": round(float(readlen.get("MajorityPctDetected")) * 100.0, 2),
                "consensusreadlength": int(readlen.get("ConsensusReadLength")),
            }

        max_readlen = max([d.get("consensusreadlength") for _, d in data.items()])
        readlens_to_plot = list(range(1, max_readlen + 1))

        headers = OrderedDict()
        headers["consensusreadlength"] = {
            "title": "Read Length (bp)",
            "description": "Predicted read length from ngsderive.",
            "scale": False,
            "format": "{:,.d}",
        }
        headers["majoritypctdetected"] = {
            "title": "Read Length: % Supporting",
            "description": "Percentage of reads which were measured at the predicted read length.",
        }
        self.general_stats_addcols(data, headers)

        samples = []
        heatdata = []

        for sample, d in data.items():
            samples.append(sample)
            _this_samples_data = {}
            for parts in d.get("evidence").split(";"):
                (k, v) = parts.split("=")
                _this_samples_data[int(k)] = int(v)

            reads = [_this_samples_data.get(this_readlen, 0) for this_readlen in readlens_to_plot]
            total_reads = sum(reads)
            reads = [r / total_reads for r in reads]
            heatdata.append(reads)

        # Config for the plot
        pconfig = {
            "id": "ngsderive_readlen_plot",
            "title": "ngsderive: Read Length",
            "xTitle": "% Evidence for Read Length",
            "yTitle": "Sample",
            "square": False,
            "legend": False,
            "datalabels": False,
        }

        self.add_section(
            name="Read length",
            anchor="ngsderive-readlen",
            description="""Predicted read length provided by ngsderive. For more information, please see
            [the documentation](https://stjudecloud.github.io/ngsderive/subcommands/readlen/).""",
            plot=heatmap.plot(heatdata, readlens_to_plot, samples, pconfig=pconfig),
        )
