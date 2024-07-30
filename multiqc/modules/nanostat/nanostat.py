import logging
from typing import List

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table
from multiqc.utils import mqc_colour

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    _KEYS_MAPPING = {
        "number_of_reads": "Number of reads",
        "number_of_bases": "Total bases",
        "number_of_bases_aligned": "Total bases aligned",
        "fraction_bases_aligned": "Fraction of bases aligned",
        "median_read_length": "Median read length",
        "mean_read_length": "Mean read length",
        "read_length_stdev": "STDEV read length",
        "n50": "Read length N50",
        "average_identity": "Average percent identity",
        "median_identity": "Median percent identity",
        "active_channels": "Active channels",
        "mean_qual": "Mean read quality",
        "median_qual": "Median read quality",
    }

    _stat_types = ("aligned", "seq summary", "fastq", "fasta", "unrecognized")

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="NanoStat",
            anchor="nanostat",
            href=["https://github.com/wdecoster/nanostat/", "https://github.com/wdecoster/nanoplot/"],
            info="Reports various statistics for long read dataset in FASTQ, BAM, or albacore sequencing summary "
            "format (supports NanoPack; NanoPlot, NanoComp).",
            extra="Programs are part of the NanoPack family for summarising results of sequencing on Oxford Nanopore "
            "methods (MinION, PromethION etc.)",
            doi="10.1093/bioinformatics/bty149",
        )

        # Find and load any NanoStat reports
        self.nanostat_data = dict()
        self.has_qscores = False
        self.has_aligned = False
        self.has_seq_summary = False
        self.has_fastq = False
        self.has_fasta = False
        for f in self.find_log_files("nanostat", filehandles=True):
            self.parse_nanostat_log(f)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        for f in self.find_log_files("nanostat/legacy", filehandles=True):
            self.parse_legacy_nanostat_log(f)

        # Filter to strip out ignored sample names
        self.nanostat_data = self.ignore_samples(self.nanostat_data)

        if len(self.nanostat_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.nanostat_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.nanostat_data, "multiqc_nanostat")

        # Stats Tables
        if self.has_aligned:
            self.nanostat_stats_table("aligned")
        if self.has_seq_summary:
            self.nanostat_stats_table("seq summary")
        if self.has_fastq:
            self.nanostat_stats_table("fastq")
        if self.has_fasta:
            self.nanostat_stats_table("fasta")

        # Quality distribution Plot
        if self.has_qscores:
            self.reads_by_quality_plot()

    def parse_nanostat_log(self, f):
        """Parse output from NanoStat

        Note: Tool can be run in two different modes, giving two variants to the output.
        To avoid overwriting keys from different modes, keys are given a suffix.
        """
        nano_stats = {}
        for line in f["f"]:
            parts = line.strip().split()
            if len(parts) == 2 and parts[0] in self._KEYS_MAPPING.keys():
                key = self._KEYS_MAPPING.get(parts[0])
                if key:
                    nano_stats[key] = float(parts[1])
            else:
                parts = line.strip().split(":")
                key = parts[0].replace("Reads ", "")
                if key.startswith(">Q"):
                    # Number of reads above Q score cutoff
                    val = int(parts[1].strip().split()[0])
                    nano_stats[key] = val
        self.save_data(f, nano_stats)

    def parse_legacy_nanostat_log(self, f):
        """Parse legacy output from NanoStat

        Note: Tool can be run in two different modes, giving two variants to the output.
        To avoid overwriting keys from different modes, keys are given a suffix.
        """
        nano_stats = {}
        for line in f["f"]:
            parts = line.strip().split(":")
            if len(parts) == 0:
                continue

            key = parts[0]

            if key in self._KEYS_MAPPING.values():
                val = float(parts[1].replace(",", ""))
                nano_stats[key] = val
            elif key.startswith(">Q"):
                # Number of reads above Q score cutoff
                val = int(parts[1].strip().split()[0])
                nano_stats[key] = val
        self.save_data(f, nano_stats)

    def save_data(self, f, nano_stats):
        """
        Normalise fields and save parsed data.

        Used for both legacy and new data formats.
        """
        if any(k.startswith(">Q") for k in nano_stats):
            self.has_qscores = True

        if "Total bases aligned" in nano_stats:
            stat_type = "aligned"
            self.has_aligned = True
        elif "Active channels" in nano_stats:
            stat_type = "seq summary"
            self.has_seq_summary = True
        elif "Mean read quality" in nano_stats:
            stat_type = "fastq"
            self.has_fastq = True
        elif "Mean read length" in nano_stats:
            stat_type = "fasta"
            self.has_fasta = True
        else:
            log.debug(f"Did not recognise NanoStat file '{f['fn']}' - skipping")
            return

        out_d = {f"{k}_{stat_type}": v for k, v in nano_stats.items()}

        # Warn if we find overlapping data for the same sample
        if f["s_name"] in self.nanostat_data:
            # Only if the same has some keys in common
            if not set(self.nanostat_data[f["s_name"]].keys()).isdisjoint(out_d.keys()):
                log.debug(f"Duplicate sample data found! Overwriting: {f['s_name']}")

        self.nanostat_data.setdefault(f["s_name"], {}).update(out_d)

        self.add_data_source(f)

    def nanostat_stats_table(self, stat_type):
        """Take the parsed stats from the Kallisto report and add it to the
        basic stats table at the top of the report"""

        headers_base = {
            "Active channels": {
                "title": "Active channels",
                "description": "Active channels",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "Median read length": {
                "title": "Median length",
                "description": "Median read length (bp)",
                "suffix": " bp",
                "format": "{:,.0f}",
                "shared_key": "nucleotides",
                "scale": "BuPu",
            },
            "Mean read length": {
                "title": "Mean length",
                "description": "Mean read length (bp)",
                "suffix": " bp",
                "scale": "Purples",
                "format": "{:,.0f}",
                "shared_key": "nucleotides",
                "hidden": True,
            },
            "Read length N50": {
                "title": "Read N50",
                "description": "Read length N50",
                "format": "{:,.0f}",
                "suffix": " bp",
                "scale": "RdPu",
            },
            "Median read quality": {
                "title": "Median Qual",
                "description": "Median read quality (Phred scale)",
                "shared_key": "phred_score",
                "scale": "RdYlGn",
            },
            "Mean read quality": {
                "title": "Mean Qual",
                "description": "Mean read quality (Phred scale)",
                "scale": "PiYG",
                "shared_key": "phred_score",
                "hidden": True,
            },
            "Median percent identity": {
                "title": "Median Identity",
                "description": "Median percent identity",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlBu",
                "shared_key": "percent_identity",
            },
            "Average percent identity": {
                "title": "Mean Identity",
                "description": "Average percent identity",
                "max": 100,
                "suffix": "%",
                "scale": "Spectral",
                "shared_key": "percent_identity",
                "hidden": True,
            },
            "Number of reads": {
                "title": f"# Reads ({config.long_read_count_prefix})",
                "description": f"Number of reads ({config.long_read_count_desc})",
                "modify": lambda x: x * config.long_read_count_multiplier,
                "shared_key": "long_read_count",
                "scale": "YlGn",
            },
            "Total bases": {
                "title": f"Total Bases ({config.base_count_prefix})",
                "description": f"Total bases ({config.base_count_desc})",
                "modify": lambda x: x * config.base_count_multiplier,
                "shared_key": "base_count",
                "scale": "BrBG",
            },
            "Total bases aligned": {
                "title": f"Aligned Bases ({config.base_count_prefix})",
                "description": f"Total bases aligned ({config.base_count_desc})",
                "modify": lambda x: x * config.base_count_multiplier,
                "shared_key": "base_count",
                "scale": "PuOr",
            },
        }

        # Add the stat_type suffix
        headers = {}
        for k in headers_base:
            key = f"{k}_{stat_type}"
            headers[key] = headers_base.get(k, dict()).copy()

        # Table config
        table_config = {
            "namespace": "NanoStat",
            "id": f"nanostat_{stat_type.replace(' ', '_')}_stats_table",
            "title": f"NanoStat {stat_type}",
        }

        # Add the report section
        self.add_section(
            name="Summary Statistics",
            anchor=f"nanostat_{stat_type.replace(' ', '_')}_stats",
            plot=table.plot(self.nanostat_data, headers, table_config),
        )

    def reads_by_quality_plot(self):
        def _get_total_reads(_data_dict):
            for _stat_type in self._stat_types:
                total_key = f"Number of reads_{_stat_type}"
                if total_key in _data_dict:
                    return _data_dict[total_key], _stat_type
            return None, None

        q_values: List[int] = []
        samples_with_same_q_values = []
        samples_with_other_q_values = []  # will be skipped
        for s_name, data_dict in self.nanostat_data.items():
            sample_q_vals = [
                int(k.strip().replace(">Q", "").split("_")[0]) for k in data_dict.keys() if k.startswith(">Q")
            ]
            if q_values and set(sample_q_vals) != set(q_values):
                samples_with_other_q_values.append(s_name)
                continue
            if not q_values:
                q_values = sample_q_vals
            samples_with_same_q_values.append(s_name)

        bar_data = {}
        q_ranges = {}
        prev_q = 0
        for q in q_values:
            if prev_q == 0:
                q_ranges[f">Q{q}"] = f"<Q{q}"
            else:
                q_ranges[f">Q{q}"] = f"Q{prev_q}-{q}"
            prev_q = q
        q_ranges["rest"] = f">Q{prev_q}"

        for s_name in samples_with_same_q_values:
            data_dict = self.nanostat_data[s_name]
            reads_total, stat_type = _get_total_reads(data_dict)
            if stat_type == "fasta":
                log.debug(f"Sample '{s_name}' has no quality metrics - excluded from quality plot")
                continue
            if s_name in bar_data and stat_type == "aligned":
                log.debug(f"Sample '{s_name}' duplicated in the quality plot - ignoring aligned data")
                continue
            elif s_name in bar_data and stat_type == "seq summary":
                log.debug(f"Sample '{s_name}' duplicated in the quality plot - overwriting with seq summary data")
            bar_data[s_name] = {}

            prev_reads = reads_total
            for k, range_name in q_ranges.items():
                if k != "rest":
                    data_key = f"{k}_{stat_type}"
                    reads_gt = data_dict[data_key]

                    bar_data[s_name][range_name] = prev_reads - reads_gt

                    if bar_data[s_name][range_name] < 0:
                        log.error(f"Error on {s_name} {range_name} {data_key}: negative number of reads")
                    prev_reads = reads_gt
                else:
                    data_key = f">Q{prev_q}_{stat_type}"
                    bar_data[s_name][range_name] = data_dict.get(data_key, 0)  # Handle 'rest' data

        cats = {}
        keys = reversed(list(q_ranges.values()))
        colours = mqc_colour.mqc_colour_scale("RdYlGn-rev", 0, len(q_ranges))
        for idx, k in enumerate(keys):
            cats[k] = {"name": "Reads " + k, "color": colours.get_colour(idx, lighten=1)}

        # Config for the plot
        pconfig = {
            "id": "nanostat_quality_dist",
            "title": "NanoStat: Reads by quality",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        content_before_plot = ""
        if samples_with_other_q_values:
            log.error(
                f"Quality score cutoffs differ between samples {samples_with_same_q_values} and {samples_with_other_q_values}, "
                f"will use values from the former: {q_values}"
            )
            content_before_plot = (
                f'<div class="alert alert-warning">Different quality score cutoffs used across different logs. '
                f'Using ones from the first met sample, and hiding samples '
                f'{", ".join("<code>" + s + "</code>" for s in samples_with_other_q_values)}.</div>'
            )

        # Add the report section
        self.add_section(
            name="Reads by quality",
            anchor="nanostat_read_qualities",
            description="Read counts categorised by read quality (Phred score).",
            content_before_plot=content_before_plot,
            helptext="""
            Sequencing machines assign each generated read a quality score using the
            [Phred scale](https://en.wikipedia.org/wiki/Phred_quality_score).
            The phred score represents the liklelyhood that a given read contains errors.
            High quality reads have a high score.
            """,
            plot=bargraph.plot(bar_data, cats, pconfig),
        )
