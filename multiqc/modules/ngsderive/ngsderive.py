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
from multiqc.plots import bargraph, heatmap

# Initialise the logger
log = logging.getLogger(__name__)

BUFFER_SIZE = 4096
COMMON_DELIMITERS = ["\t", ",", " ", "|"]
EXPECTED_HEADER_COUNT_STRANDEDNESS = 5
EXPECTED_HEADER_COUNT_INSTRUMENT = 4
EXPECTED_HEADER_COUNT_READLEN = 4


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
            info="is a forensic analysis tool useful in backwards computing information from "
            + "next-generation sequencing data. "
            + "Notably, results are provided as a 'best guess' — the tool does "
            + "not claim 100%% accuracy and results should be considered with that understanding. "
            + "Please see the <a href='https://github.com/claymcleod/ngsderive/#ngsderive'>documentation</a> "
            + "for more information. ",
        )

        self.strandedness = {}
        self.instrument = {}
        self.readlen = {}

        # parse ngsderive summary file
        for f in self.find_log_files("ngsderive/strandedness"):
            self.parse(
                self.strandedness,
                f,
                "strandedness",
                EXPECTED_HEADER_COUNT_STRANDEDNESS,
            )

        for f in self.find_log_files("ngsderive/instrument"):
            self.parse(
                self.instrument,
                f,
                "instrument",
                EXPECTED_HEADER_COUNT_INSTRUMENT,
            )

        for f in self.find_log_files("ngsderive/readlen"):
            self.parse(
                self.readlen,
                f,
                "readlen",
                EXPECTED_HEADER_COUNT_READLEN,
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

        for delim in COMMON_DELIMITERS:
            if delim in header and len(header.split(delim)) == expected_header_count:
                return {"delimiter": delim}

        log.warn("Could not easily detect delimiter for file. Trying csv.Sniffer()")
        dialect = csv.Sniffer().sniff(f.read(BUFFER_SIZE))
        f.seek(0)
        return {"dialect": dialect}

    def parse(self, sample_dict, found_file, subcommand, expected_header_count):
        contents = found_file.get("f")
        root = found_file.get("root")
        source = os.path.abspath(os.path.join(found_file["root"], found_file["fn"]))

        kwargs = self.probe_file_for_dictreader_kwargs(
            io.StringIO(contents), expected_header_count
        )

        relevant_items = []
        for row in csv.DictReader(io.StringIO(contents), **kwargs):
            if not row.get("File"):
                continue
            sample_name = self.clean_s_name(row.get("File"), root)
            if sample_name in sample_dict:
                log.debug(
                    "Duplicate sample name found! Overwriting: {}".format(sample_name)
                )

            sample_dict[sample_name] = row
            self.add_data_source(s_name=sample_name, source=source)

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
            "title": "Strandedness (Predicted)",
            "description": "Predicted strandedness from ngsderive",
        }
        headers["forward"] = {
            "title": "Strandedness (Forward %)",
            "description": "Percentage of reads which were evidence for forward-stranded.",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
            "hidden": True,
        }
        headers["reverse"] = {
            "title": "Strandedness (Reverse %)",
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
            "title": "ngsderive: strandedness",
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
            description="Predicted strandedness provided by ngsderive. "
            + 'For more information, please see <a href="https://stjudecloud.github.io/ngsderive/subcommands/strandedness/">'
            + "the relevant documentation and limitations</a>.",
            plot=bargraph.plot(bardata, headers, pconfig),
        )

    def add_instrument_data(self):
        data = {}
        for sample, instrument in self.instrument.items():
            data[sample] = {
                "instrument": "/".join(
                    sorted(instrument.get("Instrument").split(" or "))
                ),
                "confidence": instrument.get("Confidence"),
                "basis": instrument.get("Basis"),
            }

        headers = OrderedDict()
        headers["instrument"] = {
            "title": "Instrument (Predicted)",
            "description": "Predicted instrument from ngsderive",
        }
        headers["confidence"] = {
            "title": "Instrument (Confidence)",
            "description": "Level of confidence (low, medium, high) that the predicted instrument is correct.",
        }
        headers["basis"] = {
            "title": "Instrument (Basis)",
            "description": "Basis upon which the prediction was made (see documentation for more details).",
            "hidden": True,
        }
        self.general_stats_addcols(data, headers)

        samples = []
        samples_data = []
        instruments = set()

        # first pass through, start to formulate heatmap data
        for s, d in data.items():
            samples.append(s)
            data = d.get("instrument").split("/")
            samples_data.append(data)
            instruments.update(data)

        # move multiple instruments to the end if it exists
        instruments = sorted(instruments)
        if "multiple instruments" in instruments:
            instruments.remove("multiple instruments")
            instruments.append("multiple instruments")

        # one-hot encode machine for each sample
        heatdata = []
        for _sample_data in samples_data:
            heatdata.append([int(i in _sample_data) for i in instruments])

        # sort the table so that it looks nice. essentially, this is
        # done as treating the one-hot encoded vectors as a binary
        # encoded number that is then sorted.
        heatdata = zip(samples, heatdata)

        def reduce(x):
            return sum([pow(i - 1, _x) for (i, _x) in enumerate(x)])

        sorted_heatdata = sorted(heatdata, key=lambda x: reduce(x[1]))
        samples, heatdata = zip(*sorted_heatdata)

        # Config for the plot
        pconfig = {
            "id": "ngsderive_instruments_plot",
            "title": "ngsderive: instruments",
            "xTitle": "Predicted Instrument",
            "yTitle": "Sample",
            "square": False,
            "legend": False,
            "datalabels": False,
        }

        self.add_section(
            name="Instrument",
            anchor="ngsderive-instrument",
            description="Predicted instrument provided by ngsderive. "
            + 'For more information, please see <a href="https://stjudecloud.github.io/ngsderive/subcommands/instrument/">'
            + "the relevant documentation and limitations</a>.",
            plot=heatmap.plot(heatdata, instruments, samples, pconfig=pconfig),
        )

    def add_readlen_data(self):
        data = {}
        for sample, readlen in self.readlen.items():
            data[sample] = {
                "evidence": readlen.get("Evidence"),
                "majoritypctdetected": round(
                    float(readlen.get("MajorityPctDetected")) * 100.0, 2
                ),
                "consensusreadlength": int(readlen.get("ConsensusReadLength")),
            }

        max_readlen = max([d.get("consensusreadlength") for _, d in data.items()])
        readlens_to_plot = list(range(1, max_readlen + 1))

        headers = OrderedDict()
        headers["consensusreadlength"] = {
            "title": "Read Length (Predicted)",
            "description": "Predicted read length from ngsderive.",
            "scale": False,
            "format": "{:,.d}",
        }
        headers["majoritypctdetected"] = {
            "title": "Read Length (% Supporting Reads)",
            "description": "Percentage of reads which were measured at the predicted read length.",
        }
        headers["evidence"] = {
            "title": "Read Length (Evidence)",
            "description": "Evidence from readlen analysis (see documentation for more details). "
            + "Typically, this is not helpful to the user in text form.",
            "hidden": True,
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

            reads = [
                _this_samples_data.get(this_readlen, 0)
                for this_readlen in readlens_to_plot
            ]
            total_reads = sum(reads)
            reads = [r / total_reads for r in reads]
            heatdata.append(reads)

        # Config for the plot
        pconfig = {
            "id": "ngsderive_readlen_plot",
            "title": "ngsderive: readlen",
            "xTitle": "% Evidence for Read Length",
            "yTitle": "Sample",
            "square": False,
            "legend": False,
            "datalabels": False,
        }

        self.add_section(
            name="Read length",
            anchor="ngsderive-readlen",
            description="Predicted read length provided by ngsderive. "
            + 'For more information, please see <a href="https://stjudecloud.github.io/ngsderive/subcommands/readlen/">'
            + "the relevant documentation and limitations</a>.",
            plot=heatmap.plot(heatdata, readlens_to_plot, samples, pconfig=pconfig),
        )
