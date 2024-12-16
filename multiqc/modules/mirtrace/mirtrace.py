import json
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="miRTrace",
            anchor="mirtrace",
            href="https://github.com/friedlanderlab/mirtrace",
            info="Quality control for small RNA sequencing data.",
            extra="""
            miRTrace performs adapter trimming and discards the reads that fail to pass
            the QC filters. miRTrace specifically addresses sequencing quality, read length,
            sequencing depth and miRNA complexity and also identifies the presence of both
            miRNAs and undesirable sequences derived from tRNAs, rRNAs, or Illumina artifact
            sequences.
        
            miRTrace also profiles clade-specific miRNAs based on a comprehensive catalog
            of clade-specific miRNA families identified previously. With this information,
            miRTrace can detect exogenous miRNAs, which could be contamination derived,
            e.g. index mis-assignment on sample demultiplexing, or biologically derived,
            e.g. parasitic RNAs.
            """,
            doi="10.1186/s13059-018-1588-9",
        )

        # Find and load miRTrace summary statistics table
        self.summary_data = dict()
        for f in self.find_log_files("mirtrace/summary"):
            self.parse_summary(f)

        # Find and load miRTrace read length table
        self.length_data = dict()
        for f in self.find_log_files("mirtrace/length"):
            self.parse_length(f)

        # Find and load miRTrace contamination statistics summary_table
        self.contamination_data = dict()
        for f in self.find_log_files("mirtrace/contaminationbasic"):
            self.parse_contamination(f)

        # Find and load miRTrace miRNA complexity table
        self.complexity_data = dict()
        for f in self.find_log_files("mirtrace/mirnacomplexity"):
            self.parse_complexity(f)

        # Filter to strip out ignored sample names
        self.summary_data = self.ignore_samples(self.summary_data)
        self.length_data = self.ignore_samples(self.length_data)
        self.contamination_data = self.ignore_samples(self.contamination_data)
        self.complexity_data = self.ignore_samples(self.complexity_data)

        # Warning when no files are found
        if (
            max(len(self.summary_data), len(self.length_data), len(self.contamination_data), len(self.complexity_data))
            == 0
        ):
            raise ModuleNoSamplesFound

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed data to a file
        self.write_data_file(self.summary_data, "multiqc_mirtrace_summary")
        self.write_data_file(self.length_data, "multiqc_mirtrace_length")
        self.write_data_file(self.contamination_data, "multiqc_mirtrace_contamination")
        self.write_data_file(self.complexity_data, "multiqc_mirtrace_complexity")

        # Report sections
        if len(self.summary_data) > 0:
            self.add_section(name="QC Plot", anchor="mirtrace_qc", plot=self.mirtrace_qc_plot())
            self.add_section(
                name="RNA Categories", anchor="mirtrace_rna_categories", plot=self.mirtrace_rna_categories()
            )

        if len(self.length_data) > 0:
            self.add_section(
                name="Read Length Distribution", anchor="mirtrace_length", plot=self.mirtrace_length_plot()
            )

        if len(self.contamination_data) > 0:
            self.add_section(
                name="Contamination Check",
                anchor="mirtrace_contamination_check",
                plot=self.mirtrace_contamination_check(),
            )

        if len(self.complexity_data) > 0:
            self.add_section(
                name="miRNA Complexity", anchor="mirtrace_complexity", plot=self.mirtrace_complexity_plot()
            )

    # Parse a miRTrace results.json file
    def parse_summary(self, f):
        try:
            cdict = json.loads(f["f"])
        except ValueError as e:
            raise e

        if "results" in cdict.keys():
            for record in cdict["results"]:
                s_name = self.clean_s_name(record["verbosename"], f)
                parsed_data = {}
                parsed_data["filename"] = record["filename"]
                parsed_data["reads_total"] = record["stats"]["allSeqsCount"]
                parsed_data["adapter_removed_length_ok"] = record["stats"]["statsQC"][4]
                parsed_data["adapter_not_detected"] = record["stats"]["statsQC"][3]
                parsed_data["length_shorter_than_18"] = record["stats"]["statsQC"][2]
                parsed_data["low_complexity"] = record["stats"]["statsQC"][1]
                parsed_data["low_phred"] = record["stats"]["statsQC"][0]
                parsed_data["reads_mirna"] = record["stats"]["statsRNAType"][0]
                parsed_data["reads_rrna"] = record["stats"]["statsRNAType"][1]
                parsed_data["reads_trna"] = record["stats"]["statsRNAType"][2]
                parsed_data["reads_artifact"] = record["stats"]["statsRNAType"][3]
                parsed_data["reads_unknown"] = record["stats"]["statsRNAType"][4]
                if s_name in self.summary_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name)
                self.summary_data[s_name] = parsed_data
        else:
            log.debug(f"No valid data {f['fn']} in miRTrace summary")
            return None

    # Parse a miRTrace mirtrace-stats-length.tsv file
    def parse_length(self, f):
        header = []
        body = {}
        lines = f["f"].splitlines()
        for line in lines:
            s = line.split("\t")
            if len(header) == 0:
                if s[0] != "LENGTH":
                    log.debug(f"No valid data {f['fn']} for read length distribution")
                    return None
                header = s[1:]
            else:
                body[s[0]] = s[1 : len(s)]

        for record in header[0 : len(header)]:
            s_name = self.clean_s_name(record, f)
            parsed_data = {}
            idx = header[0 : len(header)].index(record)
            for length in body:
                parsed_data[length] = int(body[length][idx])
            if s_name in self.length_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.length_data[s_name] = parsed_data

    # Parse a miRTrace mirtrace-stats-contamination_basic.tsv file
    def parse_contamination(self, f):
        header = []
        body = {}
        lines = f["f"].splitlines()
        for line in lines:
            s = line.split("\t")
            if len(header) == 0:
                if s[0] != "CLADE":
                    log.debug(f"No valid data {f['fn']} for contamination check")
                    return None
                header = s[1:]
            else:
                body[s[0]] = s[1 : len(s)]

        for record in header[0 : len(header)]:
            s_name = self.clean_s_name(record, f)
            parsed_data = {}
            idx = header[0 : len(header)].index(record)
            for clade in body:
                parsed_data[clade] = int(body[clade][idx])
            if s_name in self.contamination_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.contamination_data[s_name] = parsed_data

    # Parse a miRTrace mirtrace-stats-mirna-complexity.tsv file
    def parse_complexity(self, f):
        header = []
        body = {}
        lines = f["f"].splitlines()
        for line in lines:
            s = line.split("\t")
            if len(header) == 0:
                if s[0] != "DISTINCT_MIRNA_HAIRPINS_ACCUMULATED_COUNT":
                    log.debug(f"No valid data {f['fn']} for miRNA complexity")
                    return None
                header = s[1:]
            else:
                body[s[0]] = s[1 : len(s)]

        for record in header[0 : len(header)]:
            s_name = self.clean_s_name(record, f)
            parsed_data = {}
            idx = header[0 : len(header)].index(record)
            for depth in body:
                parsed_data[depth] = int(body[depth][idx]) if body[depth][idx] else 0
            if s_name in self.complexity_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.complexity_data[s_name] = parsed_data

    # miRTrace QC Plot
    def mirtrace_qc_plot(self):
        """Generate the miRTrace QC Plot"""

        # Specify the order of the different possible categories
        keys = {
            "adapter_removed_length_ok": {"color": "#006837", "name": "Reads ≥ 18 nt after adapter removal"},
            "adapter_not_detected": {"color": "#66bd63", "name": "Reads without adapter"},
            "length_shorter_than_18": {"color": "#fdae61", "name": "Reads < 18 nt after adapter removal"},
            "low_complexity": {"color": "#d73027", "name": "Reads with low complexity"},
            "low_phred": {"color": "#a50026", "name": "Reads with low PHRED score"},
        }

        # Config for the plot
        config = {
            "id": "mirtrace_qc_plot",
            "title": "miRTrace: QC Plot",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(self.summary_data, keys, config)

    # miRTrace Read Length Distribution
    def mirtrace_length_plot(self):
        """Generate the miRTrace Read Length Distribution"""

        data = dict()
        for s_name in self.length_data:
            try:
                data[s_name] = {int(d): int(self.length_data[s_name][d]) for d in self.length_data[s_name]}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug("No valid data for read length distribution")
            return None

        config = {
            "id": "mirtrace_length_plot",
            "title": "miRTrace: Read Length Distribution",
            "ylab": "Read Count",
            "xlab": "Read Lenth (bp)",
            "ymin": 0,
            "xmin": 0,
            "x_decimals": False,
            "tt_label": "<b>Read Length (bp) {point.x}</b>: {point.y} Read Count",
            "x_bands": [
                {"from": 40, "to": 50, "color": "#ffebd1"},
                {"from": 26, "to": 40, "color": "#e2f5ff"},
                {"from": 18, "to": 26, "color": "#e5fce0"},
                {"from": 0, "to": 18, "color": "#ffffe2"},
            ],
        }

        return linegraph.plot(data, config)

    # miRTrace RNA Categories
    def mirtrace_rna_categories(self):
        """Generate the miRTrace RNA Categories"""

        # Specify the order of the different possible categories
        keys = {
            "reads_mirna": {"color": "#33a02c", "name": "miRNA"},
            "reads_rrna": {"color": "#ff7f00", "name": "rRNA"},
            "reads_trna": {"color": "#1f78b4", "name": "tRNA"},
            "reads_artifact": {"color": "#fb9a99", "name": "Artifact"},
            "reads_unknown": {"color": "#d9d9d9", "name": "Unknown"},
        }

        # Config for the plot
        config = {
            "id": "mirtrace_rna_categories_plot",
            "title": "miRTrace: RNA Categories",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(self.summary_data, keys, config)

    # miRTrace Contamination Check
    def mirtrace_contamination_check(self):
        """Generate the miRTrace Contamination Check"""

        # A library of 24 colors. Should be enough for this plot
        color_lib = [
            "#A6CEE3",
            "#1F78B4",
            "#B2DF8A",
            "#33A02C",
            "#FB9A99",
            "#E31A1C",
            "#FDBF6F",
            "#FF7F00",
            "#CAB2D6",
            "#6A3D9A",
            "#FFFF99",
            "#B15928",
            "#8DD3C7",
            "#FFFFB3",
            "#BEBADA",
            "#FB8072",
            "#80B1D3",
            "#FDB462",
            "#B3DE69",
            "#FCCDE5",
            "#D9D9D9",
            "#BC80BD",
            "#CCEBC5",
            "#FFED6F",
        ]

        idx = 0

        # Specify the order of the different possible categories
        keys = {}
        for clade in self.contamination_data[list(self.contamination_data.keys())[0]]:
            keys[clade] = {"color": color_lib[idx], "name": clade}
            if idx < 23:
                idx += 1
            else:
                idx = 0

        # Config for the plot
        config = {
            "cpswitch_c_active": False,
            "id": "mirtrace_contamination_check_plot",
            "title": "miRTrace: Contamination Check",
            "ylab": "# miRNA detected",
            "cpswitch_counts_label": "Number of detected miRNA",
        }

        return bargraph.plot(self.contamination_data, keys, config)

    # miRTrace Read Length Distribution
    def mirtrace_complexity_plot(self):
        """Generate the miRTrace miRNA Complexity Plot"""

        data = dict()
        for s_name in self.complexity_data:
            try:
                data[s_name] = {int(self.complexity_data[s_name][d]): int(d) for d in self.complexity_data[s_name]}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug("No valid data for miRNA complexity")
            return None

        config = {
            "id": "mirtrace_complexity_plot",
            "title": "miRTrace: miRNA Complexity Plot",
            "ylab": "Distinct miRNA Count",
            "xlab": "Number of Sequencing Reads",
            "ymin": 0,
            "xmin": 1,
            "x_decimals": False,
            "tt_label": "<b>Number of Sequencing Reads {point.x}</b>: {point.y} Distinct miRNA Count",
        }

        return linegraph.plot(data, config)
