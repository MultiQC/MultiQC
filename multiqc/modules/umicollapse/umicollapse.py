import logging
import re
from typing import Dict, Union

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, violin

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Sample names are extracted from log files if possible. In logs, the command
    line arguments are printed, which must have both the input and output
    file paths.

    ```
    umicollapse bam -i SRR19887568.sorted.bam -o SRR19887568.umi_dedup.sorted.bam
    Arguments	[bam, -i, SRR19887568.sorted.bam, -o, SRR19887568.umi_dedup.sorted.bam]
    ```

    `umicollapse` requires both -i and -o flags as valid file paths. Process
    substitution is not supported currently by umicollapse, but in case it is
    used for the -i flag, we fallback to the log file name.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="UMICollapse",
            anchor="umicollapse",
            href="https://github.com/Daniel-Liu-c0deb0t/UMICollapse",
            info="Algorithms for efficiently collapsing reads with Unique Molecular Identifiers",
            doi="10.7717/peerj.8275",
        )

        data_by_sample = dict()
        for f in self.find_log_files("umicollapse"):
            # Parse the log file for sample name and statistics
            data = self.parse_logs(f)
            if data:
                f["s_name"] = self._parse_s_name(f)
                if f["s_name"] in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                data_by_sample[f["s_name"]] = data
                self.add_data_source(f)  # Add the log file information to the multiqc_sources.txt
                if "version" in data:  # Add version info
                    self.add_software_version(data.pop("version"), f["s_name"])

        # Check the log files against the user supplied list of samples to ignore
        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0 and len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} UMICollapse reports")

        # Write parsed reports data to a file
        self.write_data_file(data_by_sample, "multiqc_umicollapse")

        if data_by_sample:
            self.general_stats_table(data_by_sample)
            self.bar_plot(data_by_sample)
            self.violin_plot(data_by_sample)

    def _parse_s_name(self, f) -> str:
        in_m = re.search(r"Arguments\t\[(\S+, )+-i, ([^,]+)", f["f"])
        if in_m:
            in_name = in_m.group(2)
            return self.clean_s_name(in_name, f)

        out_m = re.search(r"Arguments\t\[(\S+, )+-o, ([^,]+)", f["f"])
        if out_m:
            out_name = out_m.group(2)
            return self.clean_s_name(out_name, f)

        return f["s_name"]

    @staticmethod
    def parse_logs(f) -> Dict:
        regexes = [
            (int, "input_reads", r"Number of input reads\t(\d+)"),
            (
                int,
                "dedup_input_reads",
                r"Number of unremoved reads\t(\d+)|Number of unique reads\t(\d+)",
            ),  # BAM mode | FastQ mode
            (int, "positions_deduplicated", r"Number of unique alignment positions\t(\d+)"),
            (float, "mean_umi_per_pos", r"Average number of UMIs per alignment position\t([\d\.]+)"),
            (int, "max_umi_per_pos", r"Max number of UMIs over all alignment positions\t(\d+)"),
            (int, "dedup_output_reads", r"Number of reads after deduplicating\t(\d+)"),
        ]

        data: Dict[str, Union[str, int, float]] = {}
        # Search for values using regular expressions
        for type_, key, regex in regexes:
            comp_regex = re.compile(regex)
            re_matches = comp_regex.search(f["f"])
            if re_matches:
                data[key] = type_(re_matches.group(1) or re_matches.group(2))  # BAM mode | FastQ mode

        # Calculate a few simple supplementary stats
        try:
            data["dedup_percent_passing"] = round(((data["dedup_output_reads"] / data["dedup_input_reads"]) * 100.0), 2)
        except (KeyError, ZeroDivisionError):
            pass
        try:
            data["removed_reads"] = data["input_reads"] - data["dedup_input_reads"]
        except KeyError:
            pass
        try:
            data["dedup_removed_reads"] = data["dedup_input_reads"] - data["dedup_output_reads"]
        except KeyError:
            pass

        return data

    def general_stats_table(self, data_by_sample):
        """
        Take the parsed stats from the umicollapse report and add it to the
        basic stats table at the top of the report
        """
        headers = {
            "dedup_output_reads": {
                "title": "Unique Reads",
                "description": f"Reads remaining after deduplication ({config.read_count_desc})",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "scale": "PuRd",
            },
            "dedup_percent_passing": {
                "title": "Pass Dedup",
                "description": "% processed reads that passed deduplication",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
        }
        self.general_stats_addcols(data_by_sample, headers)

    def bar_plot(self, data_by_sample):
        """Generate a plot of the read deduplication rates for the main report"""

        # Specify the order of the different possible categories
        keys = {
            "dedup_output_reads": {"color": "#7fc97f", "name": "Reads remaining"},
            "dedup_removed_reads": {"color": "#fdc086", "name": "Reads removed"},
        }

        return self.add_section(
            name="Deduplicated Reads",
            anchor="umicollapse-dedup-plot",
            plot=bargraph.plot(
                data_by_sample,
                keys,
                {
                    "id": "umicollapse_deduplication_barplot",
                    "title": "UMI-tools: Deduplication Counts",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of Reads",
                },
            ),
        )

    def violin_plot(self, data_by_sample):
        """Generate a violin of UMI stats for the main report"""

        headers = {
            "positions_deduplicated": {
                "title": "Positions Dedup",
                "description": "Total number of positions deduplicated",
                "min": 0,
                "format": "{:,d}",
                "scale": "Greens",
            },
            "mean_umi_per_pos": {
                "title": "Mean #UMI",
                "description": "Mean number of unique UMIs per position",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Reds",
            },
            "max_umi_per_pos": {
                "title": "Max #UMI",
                "description": "Max number of unique UMIs per position",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Oranges",
            },
            "removed_reads": {
                "title": "Discarded reads",
                "description": "Discarded unmapped, unpaired, or chimeric reads",
                "min": 0,
                "format": "{:,d}",
                "scale": "Blues",
            },
        }

        # add a section with a violin plot of UMI stats to the report
        self.add_section(
            name="UMI Stats",
            anchor="umicollapse-umi-stats",
            description="Statistics from running `umicollapse`",
            helptext="""
            - **Positions Dedup**: Total number of positions deduplicated
            - **Mean #UMI**: Mean number of unique UMIs per position
            - **Max #UMI**: Max number of unique UMIs per position
            - **Discarded reads**: Discarded unmapped, unpaired, or chimeric reads
            """,
            plot=violin.plot(
                data_by_sample,
                headers,
                {
                    "id": "umicollapse_stats_violin",
                    "title": "UMICollapse: UMI stats",
                },
            ),
        )
