import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

log = logging.getLogger(__name__)

VERSION_REGEX = r"# Version: v([\d\.]+)"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Illumina InterOp Statistics",
            anchor="interop",
            href="http://illumina.github.io/interop/index.html",
            info=" - a set of common routines used for reading and writing InterOp metric files.",
            # No publication / DOI // doi=
        )

        log = logging.getLogger(__name__)

        # Parse data
        self.runSummary = {}
        self.indexSummary = {}
        for f in self.find_log_files("interop/summary", filehandles=True):
            parsed_data, version = self.parse_summary_csv(f["f"])
            if max(len(parsed_data["summary"]), len(parsed_data["details"])) > 0:
                self.runSummary[f["s_name"]] = parsed_data
                self.add_data_source(f)

            if version is not None:
                self.add_software_version(version, f["s_name"])

        for f in self.find_log_files("interop/index-summary", filehandles=True):
            parsed_data, version = self.parse_index_summary_csv(f["f"])
            if max(len(parsed_data["summary"]), len(parsed_data["details"])) > 0:
                self.indexSummary[f["s_name"]] = parsed_data
                self.add_data_source(f)

            if version is not None:
                self.add_software_version(version, f["s_name"])

        self.runSummary = self.ignore_samples(self.runSummary)
        self.indexSummary = self.ignore_samples(self.indexSummary)

        # No samples
        num_samples = max(len(self.runSummary), len(self.indexSummary))
        if num_samples == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {num_samples} reports")

        # Create report Sections
        if len(self.runSummary) > 0:
            # Write data to file
            self.write_data_file(self.runSummary, "interop_runsummary")

            self.add_section(
                name="Read Metrics Summary",
                anchor="interop-runmetrics-summary",
                description="Summary statistics for Total read count from each run.",
                plot=self.run_metrics_summary_table(self.runSummary),
            )
            self.add_section(
                name="Read Metrics per Lane",
                anchor="interop-runmetrics-details",
                plot=self.run_metrics_details_table(self.runSummary),
            )

        if len(self.indexSummary) > 0:
            # Write data to file
            self.write_data_file(self.indexSummary, "interop_indexsummary")

            self.add_section(
                name="Indexing QC Metrics summary",
                anchor="interop-indexmetrics-summary",
                description="Summary metrics about each lane",
                plot=self.index_metrics_summary_table(self.indexSummary),
            )
            self.add_section(
                name="Indexing QC Metrics details",
                anchor="interop-indexmetrics-details",
                description=" Detail Metrics about each lane",
                plot=self.index_metrics_details_table(self.indexSummary),
            )

    @staticmethod
    def parse_summary_csv(f):
        metrics = {"summary": {}, "details": {}}
        header = []
        section = None
        read = None
        version = None
        for line in f:
            line = line.strip()
            # assume fixed file format

            # find version
            if line.startswith("# Version"):
                match = re.search(VERSION_REGEX, line)
                if match:
                    version = match.group(1)

            # find summary header
            if line.startswith("Level,Yield,Projected Yield,Aligned,Error Rate,Intensity C1,%>=Q30"):
                # set section to summary
                section = "summary"
                header = line.split(",")
            elif section == "summary":
                data = line.split(",")
                # process summary
                metrics["summary"][data[0]] = {}
                for idx in range(1, len(data)):
                    try:
                        metrics["summary"][data[0]][header[idx]] = int(data[idx])
                    except ValueError:
                        try:
                            metrics["summary"][data[0]][header[idx]] = float(data[idx])
                        except ValueError:
                            metrics["summary"][data[0]][header[idx]] = data[idx]
                if line.startswith("Total"):
                    section = None
            elif line.startswith("Read") and (section is None or section == "details"):
                # set section to details
                section = "details"
                read = line
            elif line.startswith("Lane,Surface,Tiles,Density,Cluster") and section == "details":
                # get details header
                header = line.split(",")
            elif section == "details":
                if line.startswith("Extracted: "):
                    section = "finish"
                    continue
                data = line.split(",")
                # process summary
                linedata = {}
                # Check if "surface" is total (-) else skip
                if data[1] == "-":
                    metrics["details"][f"Lane {data[0]} - {read}"] = {}
                else:
                    continue
                for idx in range(2, len(data)):
                    try:
                        if header[idx] == "Phas/Prephas":
                            val = data[idx].split("/")
                            linedata["Phased"] = val[0]
                            linedata["Prephased"] = val[1]
                        elif header[idx] == "Cycles Error":
                            vals = data[idx].split(" - ")
                            linedata[header[idx]] = max(int(v) for v in vals)
                        else:
                            try:
                                linedata[header[idx]] = int(data[idx])
                            except ValueError:
                                linedata[header[idx]] = float(data[idx])
                    except ValueError:
                        val = re.sub(pattern=r"\+/-.*", repl="", string=data[idx]).strip()
                        try:
                            linedata[header[idx]] = int(val)
                        except ValueError:
                            try:
                                linedata[header[idx]] = float(val)
                            except ValueError:
                                linedata[header[idx]] = val
                metrics["details"][f"Lane {data[0]} - {read}"] = linedata

        return metrics, version

    @staticmethod
    def parse_index_summary_csv(f):
        metrics = {"summary": {}, "details": {}}
        header = []
        section = None
        lane = None

        summary = {}
        details = {}
        version = None
        for line in f:
            line = line.strip()
            # assume fixed file format
            # find version
            if line.startswith("# Version"):
                match = re.search(VERSION_REGEX, line)
                if match:
                    version = match.group(1)

            if line.startswith("Lane"):
                # set lane
                lane = line
                summary[lane] = {}
                continue
            if line.startswith("Total Reads,PF Reads,% Read Identified (PF),CV,Min,Max"):
                header = line.split(",")
                section = "summary"
                continue
            if line.startswith("Index Number,Sample Id,Project,Index 1 (I7),Index 2 (I5),% Read Identified (PF)"):
                header = line.split(",")
                section = "details"
                continue
            if section == "summary":
                data = line.split(",")
                for idx in range(0, len(data)):
                    summary[lane][header[idx]] = data[idx]
                continue
            if section == "details":
                data = line.split(",")
                details[f"{data[1]} - {lane}"] = {}
                for idx in range(2, len(data)):
                    details[f"{data[1]} - {lane}"][header[idx]] = data[idx]
                continue

        metrics["summary"] = summary
        metrics["details"] = details

        return metrics, version

    @staticmethod
    def run_metrics_summary_table(data):
        headers = {
            "Yield": {
                "rid": "summary_Yield",
                "title": "Gbp Yield",  # Numbers are rounded up to Gbp, so no point in using multiplier for smaller numbers as will be all zeroes
                "description": 'The number of bases sequenced (Gbp base pairs over all "usable cycles")',
                "scale": "PuOr",
                "min": 0,
                "format": "{:,.2f}",
            },
            "Aligned": {
                "rid": "summary_Aligned",
                "title": "Aligned (%)",
                "description": "The percentage of the sample that aligned to the PhiX genome",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "PiYG",
            },
            "Error Rate": {
                "title": "Error Rate (%)",
                "description": "",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
            },
            "Intensity C1": {
                "rid": "summary_Intensity_C1",
                "title": "Intensity Cycle 1",
                "description": "The intensity statistic at cycle 1.",
                "min": 0,
                "scale": "PuOr",
                "format": "{:,d}",
            },
            "%>=Q30": {
                "rid": "summary_Q30",
                "title": "% >= Q30",
                "description": "Percentage of reads with quality phred score of 30 or above",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
            },
        }
        table_config = {
            "namespace": "interop",
            "id": "interop-runmetrics-summary-table",
            "table_title": "Read metrics summary",
            "col1_header": "Run - Read",
        }

        tdata = {}
        for s_name in data:
            for key in data[s_name]["summary"]:
                tdata[f"{s_name} - {key}"] = data[s_name]["summary"][key]

        return table.plot(tdata, headers, table_config)

    @staticmethod
    def run_metrics_details_table(data):
        headers = {
            "Surface": {"title": "Surface", "description": ""},
            "Tiles": {"title": "Tiles", "description": "The number of tiles per lane.", "hidden": True},
            "Density": {
                "title": "Density",
                "description": "The density of clusters (in thousands per mm2) detected by image analysis, +/- 1 standard deviation.",
                "hidden": True,
            },
            "Cluster PF": {
                "title": "Cluster PF (%)",
                "description": "The percentage of clusters passing filtering, +/- 1 standard deviation.",
                "suffix": "%",
            },
            "% Occupied": {
                "title": "Occupied (%)",
                "description": "The percentage of nanowells occupied by clusters, +/- 1 standard deviation.",
                "suffix": "%",
            },
            "Phased": {
                "title": "Phased (%)",
                "description": "The value used by RTA for the percentage of molecules in a cluster for which sequencing falls behind (phasing) or jumps ahead (prephasing) the current cycle within a read.",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
            },
            "Prephased": {
                "title": "Prephased (%)",
                "description": "The value used by RTA for the percentage of molecules in a cluster for which sequencing falls behind (phasing) or jumps ahead (prephasing) the current cycle within a read.",
                "format": "{:.,2f}",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
            },
            "Reads": {
                "title": "M Reads",
                "description": "The number of clusters (millions)",
            },
            "Reads PF": {
                "title": "M PF Reads",
                "description": "The number of passing filter clusters (millions)",
            },
            "Cycles Error": {
                "title": "Cycles Error",
                "description": "The number of cycles that have been error-rated using PhiX, starting at cycle 1.",
                "format": "{:.,0f}",
                "scale": "OrRd",
            },
            "Yield": {
                "title": "Gbp Yield",
                "description": "The number of bases sequenced which passed filter (Gbp base pairs)",
                "scale": "PuOr",
                "min": 0,
                "format": "{:,.2f}",
            },
            "Aligned": {
                "title": "Aligned (%)",
                "description": "The percentage that aligned to the PhiX genome.",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "PiYG",
            },
            "Error": {
                "title": "Error Rate (%)",
                "description": "The calculated error rate, as determined by the PhiX alignment.",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
            },
            "Error (35)": {
                "title": "Error Rate 35 Cycles (%)",
                "description": "The calculated error rate for cycles 1-35.",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "hidden": True,
            },
            "Error (75)": {
                "title": "Error Rate 75 Cycles (%)",
                "description": "The calculated error rate for cycles 1-75.",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "hidden": True,
            },
            "Error (100)": {
                "title": "Error Rate 100 Cycles (%)",
                "description": "The calculated error rate for cycles 1-100.",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "hidden": True,
            },
            "Intensity C1": {
                "title": "Intensity Cycle 1",
                "description": "The intensity statistic at cycle 1.",
                "min": 0,
                "scale": "PuOr",
                "format": "{:,d}",
            },
            "%>=Q30": {
                "title": "%>=Q30",
                "description": "The percentage of bases with a quality score of 30 or higher, respectively.",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
            },
        }
        table_config = {
            "namespace": "interop",
            "id": "interop-runmetrics-detail-table",
            "table_title": "Sequencing Lane Statistics",
            "col1_header": "Run - Lane - Read",
        }

        tdata = {}
        for s_name in data:
            for key in data[s_name]["details"]:
                tdata[f"{s_name} - {key}"] = data[s_name]["details"][key]

        return table.plot(tdata, headers, table_config)

    @staticmethod
    def index_metrics_summary_table(data):
        headers = {
            "Total Reads": {
                "rid": "interop_reads_total",
                "title": f"{config.read_count_prefix} Reads",
                "description": f"The total number of reads for this lane ({config.read_count_desc})",
                "modify": lambda x: float(x) * config.read_count_multiplier,
                "format": "{:,.2f}",
                "shared_key": "read_count",
            },
            "PF Reads": {
                "title": f"{config.read_count_prefix} PF Reads",
                "description": "The total number of passing filter reads for this lane ({})".format(
                    config.read_count_desc
                ),
                "modify": lambda x: float(x) * config.read_count_multiplier,
                "format": "{:,.2f}",
                "shared_key": "read_count",
            },
            "% Read Identified (PF)": {
                "rid": "summary_reads_identified_pf",
                "title": "% Reads Identified (PF)",
                "description": "The total fraction of passing filter reads assigned to an index.",
                "suffix": "%",
            },
            "CV": {
                "title": "CV",
                "description": "The coefficient of variation for the number of counts across all indexes.",
                "format": "{:.,2f}",
            },
            "Min": {"title": "Min", "description": "The lowest representation for any index."},
            "Max": {"title": "Max", "description": "The highest representation for any index."},
        }
        table_config = {
            "namespace": "interop",
            "id": "interop-indexmetrics-summary-table",
            "table_title": "Index Read Statistics Summary",
            "col1_header": "Run - Lane",
        }

        tdata = {}
        for s_name in data:
            for key in data[s_name]["summary"]:
                tdata[f"{s_name} - {key}"] = data[s_name]["summary"][key]

        return table.plot(tdata, headers, table_config)

    @staticmethod
    def index_metrics_details_table(data):
        headers = {
            "% Read Identified (PF)": {
                "title": "% Reads Identified (PF)",
                "description": "The number of reads (only includes Passing Filter reads) mapped to this index.",
                "suffix": "%",
            },
            "Index 1 (I7)": {"title": "Index 1 (I7)", "description": "The sequence for the first Index Read."},
            "Index 2 (I5)": {"title": "Index 2 (I5)", "description": "The sequence for the second Index Read."},
        }

        table_config = {
            "namespace": "interop",
            "id": "interop-indexmetrics-details-table",
            "table_title": "Index Read Statistics Details",
            "col1_header": "Run - Sample - Lane",
        }

        tdata = {}
        for s_name in data:
            for key in data[s_name]["details"]:
                tdata[f"{s_name} - {key}"] = data[s_name]["details"][key]

        return table.plot(tdata, headers, table_config)
