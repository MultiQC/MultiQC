import csv
import io
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ngsderive",
            anchor="ngsderive",
            href="https://github.com/stjudecloud/ngsderive",
            info="Forensic tool for by backwards computing library information in sequencing data",
            extra="""
            Results are provided as a 'best guess' — the tool does not claim 100% accuracy and results 
            should be considered with that understanding. Please see the documentation for more information.
            """,
            # Can't find a DOI // doi=
        )

        self.strandedness = {}
        self.instrument = {}
        self.readlen = {}
        self.encoding = {}
        self.junctions = {}

        # parse ngsderive summary file
        expected_header_count_strandedness = 5
        expected_header_count_instrument = 4
        expected_header_count_readlen = 4
        expected_header_count_encoding = 3
        expected_header_count_junctions = 9

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

        for f in self.find_log_files("ngsderive/encoding"):
            self.parse(
                self.encoding,
                f,
                "encoding",
                expected_header_count_encoding,
            )

        for f in self.find_log_files("ngsderive/junction_annotation"):
            self.parse(
                self.junctions,
                f,
                "junctions",
                expected_header_count_junctions,
            )

        self.strandedness = self.ignore_samples(self.strandedness)
        self.instrument = self.ignore_samples(self.instrument)
        self.readlen = self.ignore_samples(self.readlen)
        self.encoding = self.ignore_samples(self.encoding)
        self.junctions = self.ignore_samples(self.junctions)

        num_results_found = max(
            [len(d) for d in [self.strandedness, self.instrument, self.readlen, self.encoding, self.junctions]]
        )
        if num_results_found == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {num_results_found} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        if self.strandedness:
            self.add_strandedness_data()

        if self.instrument:
            self.add_instrument_data()

        if self.readlen:
            self.add_readlen_data()

        if self.encoding:
            self.add_encoding_data()

        if self.junctions:
            self.add_junctions_data()

    @staticmethod
    def probe_file_for_dictreader_kwargs(f, expected_header_count):
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

        log.warning("Could not easily detect delimiter for file. Trying csv.Sniffer()")
        dialect = csv.Sniffer().sniff(f.read(4096))  # 4096 is the buffer size
        f.seek(0)
        return {"dialect": dialect}

    def parse(self, sample_dict, found_file, subcommand, expected_header_count):
        kwargs = self.probe_file_for_dictreader_kwargs(io.StringIO(found_file["f"]), expected_header_count)
        for row in csv.DictReader(io.StringIO(found_file["f"]), **kwargs):
            if not row.get("File"):
                continue
            sample_name = self.clean_s_name(row.get("File"), found_file)
            if sample_name in sample_dict:
                log.debug(f"Duplicate sample name found! Overwriting: {sample_name}")

            sample_dict[sample_name] = row
            self.add_data_source(f=found_file, s_name=sample_name)

    def add_strandedness_data(self):
        # Write data to file
        self.write_data_file(self.strandedness, "ngsderive_strandedness")

        data = {}
        for sample, strandedness in self.strandedness.items():
            data[sample] = {
                "predicted": strandedness.get("Predicted"),
                "forward": round(float(strandedness.get("ForwardPct")) * 100.0, 2),
                "reverse": round(float(strandedness.get("ReversePct")) * 100.0, 2),
            }

        bardata = {}
        sorted_data = sorted(data.items(), key=lambda x: x[1].get("forward"))
        for k, v in sorted_data:
            bardata[k] = v

        headers = {
            "predicted": {
                "title": "Strandedness",
                "description": "Predicted strandedness from ngsderive",
                "scale": False,
            }
        }
        self.general_stats_addcols(data, headers)

        # Config for the plot
        pconfig = {
            "id": "ngsderive_strandedness_plot",
            "title": "ngsderive: Strandedness",
            "ylab": "% Read Evidence",
            "ymin": 0,
            "ymax": 100,
            "ysuffix": "%",
            "cpswitch": False,
        }

        self.add_section(
            name="Strandedness",
            anchor="ngsderive-strandedness",
            description="""Predicted strandedness provided by ngsderive. For more information, please see
            [the documentation](https://stjudecloud.github.io/ngsderive/subcommands/strandedness/).""",
            plot=bargraph.plot(bardata, ["forward", "reverse"], pconfig),
        )

    def add_instrument_data(self):
        # Write data to file
        self.write_data_file(self.instrument, "ngsderive_instrument")

        bgcols = {
            "low confidence": "#f8d7da",
            "medium confidence": "#fff3cd",
            "high confidence": "#d1e7dd",
        }
        cond_formatting_rules = {
            "pass": [{"s_eq": "high confidence"}],
            "warn": [{"s_eq": "medium confidence"}],
            "fail": [{"s_eq": "low confidence"}],
        }

        general_data = {}
        for sample, instrument_data in self.instrument.items():
            general_data[sample] = {
                "instrument": " / ".join(sorted(instrument_data.get("Instrument").split(" or "))),
                "confidence": instrument_data.get("Confidence"),
            }

        general_headers = {
            "instrument": {
                "title": "Predicted Instrument",
                "description": "Predicted instrument from ngsderive",
            },
            "confidence": {
                "title": "Instrument: Confidence",
                "description": "Level of confidence (low, medium, high) that the predicted instrument is correct.",
                "bgcols": bgcols,
                "cond_formatting_rules": cond_formatting_rules,
                "hidden": True,
            },
        }
        self.general_stats_addcols(general_data, general_headers)

        instruments = set()

        for d in general_data.values():
            instruments.update(d.get("instrument").split(" / "))

        # move multiple instruments to the end if it exists
        instruments = sorted(instruments)
        if "multiple instruments" in instruments:
            instruments.remove("multiple instruments")
            instruments.append("multiple instruments")

        headers = {}
        for instrument in instruments:
            headers[instrument] = {
                "title": instrument,
                "description": f"Predicted {instrument} from ngsderive",
                "bgcols": bgcols,
                "cond_formatting_rules": cond_formatting_rules,
            }
        headers["basis"] = {
            "title": "Instrument: Basis",
            "description": "Basis upon which the prediction was made.",
        }

        table_data = {}
        for sample, instrument_data in self.instrument.items():
            table_data[sample] = {}
            for instrument in instrument_data.get("Instrument").split(" or "):
                table_data[sample][instrument] = instrument_data.get("Confidence")
            table_data[sample]["basis"] = instrument_data.get("Basis")

        # Config for the plot
        config = {
            "id": "ngsderive_instruments_plot",
            "title": "ngsderive: Instruments",
            "scale": False,
        }

        self.add_section(
            name="Instrument",
            anchor="ngsderive-instrument",
            description="""Predicted instrument provided by ngsderive. For more information, please see
            [the documentation](https://stjudecloud.github.io/ngsderive/subcommands/instrument/).""",
            plot=table.plot(table_data, headers, config),
        )

    def add_readlen_data(self):
        # Write data to file
        self.write_data_file(self.readlen, "ngsderive_readlen")

        data = {}
        for sample, readlen in self.readlen.items():
            data[sample] = {
                "evidence": readlen.get("Evidence"),
                "majoritypctdetected": round(float(readlen.get("MajorityPctDetected").strip("%")) * 100.0, 2),
                "consensusreadlength": int(readlen.get("ConsensusReadLength")),
            }

        headers = {
            "consensusreadlength": {
                "title": "Read Length (bp)",
                "description": "Predicted read length from ngsderive.",
                "format": "{:,d}",
            },
            "majoritypctdetected": {
                "title": "Read Length: % Supporting",
                "description": "Percentage of reads which were measured at the predicted read length.",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "hidden": True,
            },
        }
        self.general_stats_addcols(data, headers)

        linedata = [{}, {}]

        for sample, d in data.items():
            # Build dict of count data
            count_data = {}
            for parts in d.get("evidence").split(";"):
                (k, v) = parts.split("=")
                count_data[int(k)] = int(v)
            linedata[1][sample] = count_data

            # Build dict of percentage data
            total_reads = sum(count_data.values())
            linedata[0][sample] = {readlen: (count / total_reads) * 100.0 for readlen, count in count_data.items()}

        # Config for the plot
        pconfig = {
            "id": "ngsderive_readlen_plot",
            "title": "ngsderive: Read Length",
            "xlab": "Read Length",
            "ylab": "% Evidence for Read Length",
            "xsuffix": " bp",
            "data_labels": [
                {
                    "name": "Percentages",
                    "ylab": "% Evidence for Read Length",
                    "ymax": 100,
                    "y_decimals": 2,
                },
                {"name": "Counts", "ylab": "Number of reads"},
            ],
        }

        self.add_section(
            name="Read length",
            anchor="ngsderive-readlen",
            description="""Predicted read length provided by ngsderive. For more information, please see
            [the documentation](https://stjudecloud.github.io/ngsderive/subcommands/readlen/).""",
            plot=linegraph.plot(linedata, pconfig),
        )

    def add_encoding_data(self):
        # Write data to file
        self.write_data_file(self.encoding, "ngsderive_encoding")

        general_data = {}
        for sample, encoding_data in self.encoding.items():
            general_data[sample] = {
                "probable_encoding": encoding_data.get("ProbableEncoding"),
                "evidence": encoding_data.get("Evidence"),
            }

        general_headers = {
            "probable_encoding": {
                "title": "Probable Encoding",
                "description": "Predicted PHRED score encoding from ngsderive",
            },
            "evidence": {
                "title": "Encoding: Evidence",
                "description": "Observed ASCII value ranges in PHRED score encoding",
                "hidden": True,
            },
        }
        self.general_stats_addcols(general_data, general_headers)

    def add_junctions_data(self):
        # Write data to file
        self.write_data_file(self.junctions, "ngsderive_junctions")

        data = {}
        for sample, junctions_data in self.junctions.items():
            data[sample] = {
                "total_junctions": junctions_data.get("total_junctions"),
                "known_junctions": junctions_data.get("known_junctions"),
                "partial_novel_junctions": junctions_data.get("partial_novel_junctions"),
                "novel_junctions": junctions_data.get("complete_novel_junctions"),
                "total_splice_events": junctions_data.get("total_splice_events"),
                "known_spliced_reads": junctions_data.get("known_spliced_reads"),
                "partial_novel_spliced_reads": junctions_data.get("partial_novel_spliced_reads"),
                "novel_spliced_reads": junctions_data.get("complete_novel_spliced_reads"),
            }

        headers = {
            "total_junctions": {
                "title": "Total junctions",
                "description": "Total number of junctions found by ngsderive",
                "hidden": True,
                "format": "{:,d}",
            },
            "known_junctions": {
                "title": "Known junctions",
                "description": "Number of annotated junctions found by ngsderive",
                "format": "{:,d}",
            },
            "partial_novel_junctions": {
                "title": "Partially novel junctions",
                "description": "Number of partially annotated junctions found by ngsderive",
                "format": "{:,d}",
            },
            "novel_junctions": {
                "title": "Novel junctions",
                "description": "Number of completely novel junctions found by ngsderive",
                "format": "{:,d}",
            },
            "total_splice_events": {
                "title": "Total splice events",
                "description": "Total number of spliced reads found by ngsderive",
                "hidden": True,
                "format": "{:,d}",
            },
            "known_spliced_reads": {
                "title": "Annotated spliced reads",
                "description": "Number of annotated spliced reads found by ngsderive",
                "hidden": True,
                "format": "{:,d}",
            },
            "partial_novel_spliced_reads": {
                "title": "Partially annotated spliced reads",
                "description": "Number of partially annotated spliced reads found by ngsderive",
                "hidden": True,
                "format": "{:,d}",
            },
            "novel_spliced_reads": {
                "title": "Novel spliced reads",
                "description": "Number of completely un-annotated spliced reads found by ngsderive",
                "hidden": True,
                "format": "{:,d}",
            },
        }
        self.general_stats_addcols(data, headers)

        bardata = {}
        sorted_junction_data = sorted(data.items(), key=lambda x: int(x[1].get("total_junctions")), reverse=True)
        for k, v in sorted_junction_data:
            bardata[k] = {
                "known_junctions": v["known_junctions"],
                "partial_novel_junctions": v["partial_novel_junctions"],
                "novel_junctions": v["novel_junctions"],
                "known_spliced_reads": v["known_spliced_reads"],
                "partial_novel_spliced_reads": v["partial_novel_spliced_reads"],
                "novel_spliced_reads": v["novel_spliced_reads"],
            }

        # Config for the plot
        pconfig = {
            "id": "ngsderive_junctions_plot",
            "title": "ngsderive: Junction Annotation",
            "cpswitch_counts_label": "Number",
            "y_decimals": False,
            "ylab": "Number of junctions",
            "data_labels": [
                {"name": "Junctions", "ylab": "Number of junctions"},
                {"name": "Spliced Reads", "ylab": "Number of spliced reads"},
            ],
        }

        cats = [
            {
                "known_junctions": {"name": "Known junctions"},
                "partial_novel_junctions": {"name": "Partially novel junctions"},
                "novel_junctions": {"name": "Novel junctions"},
            },
            {
                "known_spliced_reads": {"name": "Annotated spliced reads"},
                "partial_novel_spliced_reads": {"name": "Partially annotated spliced reads"},
                "novel_spliced_reads": {"name": "Novel spliced reads"},
            },
        ]
        self.add_section(
            name="Junction Annotations",
            anchor="ngsderive-junctions",
            description="""Junction annotations provided by ngsderive. For more information, please see
            [the documentation](https://stjudecloud.github.io/ngsderive/subcommands/junction_annotation/).""",
            plot=bargraph.plot([bardata, bardata], cats, pconfig),
        )
