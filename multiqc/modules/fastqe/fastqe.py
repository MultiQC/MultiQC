# https://github.com/fastqe/fastqe
# https://github.com/fastqe/fastqe/issues/11
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
from multiqc.plots.table_object import TableConfig

log = logging.getLogger(__name__)

# Emoji to Phred quality score mapping, based on FastQE default emoji encoding.
# See: https://github.com/fastqe/fastqe/blob/master/fastqe/map_data.py
EMOJI_TO_PHRED: dict[str, int] = {
    "🚫": 0,
    "❌": 1,
    "👺": 2,
    "💔": 3,
    "🙅": 4,
    "👾": 5,
    "🙈": 6,
    "🙊": 7,
    "🙉": 8,
    "🤡": 9,
    "👹": 10,
    "😡": 11,
    "😱": 12,
    "😰": 13,
    "😅": 14,
    "😭": 15,
    "😒": 16,
    "🤢": 17,
    "😩": 18,
    "😧": 19,
    "😿": 20,
    "😾": 21,
    "😐": 22,
    "😶": 23,
    "😙": 24,
    "😊": 25,
    "😃": 26,
    "😘": 27,
    "😁": 28,
    "😗": 29,
    "😚": 30,
    "😇": 31,
    "😄": 32,
    "😎": 33,
    "😉": 34,
    "😝": 35,
    "😛": 36,
    "😋": 37,
    "😌": 38,
    "😜": 39,
    "😆": 40,
    "🤗": 41,
}


class MultiqcModule(BaseMultiqcModule):
    """
    FastQE uses emoji to represent FASTQ sequence quality scores, providing a
    fun and visually intuitive way to assess sequencing data quality.

    The module parses the tab-separated output from FastQE and displays the
    emoji quality strings for each sample along with summary statistics.

    **Note** — MultiQC parses the standard output from FastQE. You must capture
    FastQE stdout to a file when running, for example:

    ```bash
    fastqe input.fastq > fastqe_output.txt
    ```
    """

    def __init__(self):
        super().__init__(
            name="FastQE",
            anchor="fastqe",
            href="https://github.com/fastqe/fastqe",
            info="Uses emoji to represent FASTQ sequence quality scores.",
        )

        fastqe_data: dict[str, dict[str, str]] = {}
        for f in self.find_log_files("fastqe", filehandles=True):
            parsed = self._parse_fastqe_log(f)
            if parsed:
                for s_name, sample_data in parsed.items():
                    s_name = self.clean_s_name(s_name, f)
                    if s_name in fastqe_data:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    fastqe_data[s_name] = sample_data
                self.add_data_source(f)

        fastqe_data = self.ignore_samples(fastqe_data)

        if len(fastqe_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(fastqe_data)} reports")

        self.add_software_version(None)

        self.add_section(
            name="Quality Scores",
            anchor="fastqe-quality",
            description="Per-base quality scores represented as emoji. "
            "Each emoji maps to a Phred quality score range.",
            plot=self._fastqe_emoji_table(fastqe_data),
        )

        numeric_data = self._compute_numeric_stats(fastqe_data)
        if numeric_data:
            headers: dict = {
                "avg_quality": {
                    "title": "Avg Quality",
                    "description": "Average Phred quality score across all bases",
                    "min": 0,
                    "max": 41,
                    "scale": "RdYlGn",
                    "format": "{:,.1f}",
                },
                "min_quality": {
                    "title": "Min Quality",
                    "description": "Minimum Phred quality score",
                    "min": 0,
                    "max": 41,
                    "scale": "RdYlGn",
                },
                "max_quality": {
                    "title": "Max Quality",
                    "description": "Maximum Phred quality score",
                    "min": 0,
                    "max": 41,
                    "scale": "RdYlGn",
                },
                "num_bases": {
                    "title": "# Bases",
                    "description": "Number of base positions",
                    "format": "{:,d}",
                    "scale": "Blues",
                },
            }
            self.general_stats_addcols(numeric_data, headers)

            self.add_section(
                name="Quality Statistics",
                anchor="fastqe-stats",
                description="Numeric quality score statistics derived from emoji encodings.",
                plot=table.plot(
                    numeric_data,
                    headers,
                    TableConfig(
                        id="fastqe_stats_table",
                        title="FastQE: Quality Statistics",
                    ),
                ),
            )

        self.write_data_file(fastqe_data, "multiqc_fastqe")

    def _parse_fastqe_log(self, f) -> dict[str, dict[str, str]] | None:
        """Parse FastQE TSV: Filename\\tStatistic\\tQualities -> {sample: {stat: emoji}}"""
        data: dict[str, dict[str, str]] = {}

        for line in f["f"]:
            line = line.strip()
            if not line:
                continue

            if line.startswith("Filename") or line.startswith("Sample Name"):
                continue

            fields = line.split("\t")
            if len(fields) >= 3:
                sample_name = fields[0].strip()
                qual_type = fields[1].strip()
                emoji_quals = fields[2].strip()
                if sample_name and qual_type and emoji_quals:
                    if sample_name not in data:
                        data[sample_name] = {}
                    data[sample_name][qual_type] = emoji_quals

        return data if data else None

    @staticmethod
    def _fastqe_emoji_table(fastqe_data: dict[str, dict[str, str]]) -> str:
        """Generate an HTML table showing emoji quality strings per sample."""
        rows = []
        for s_name, data in sorted(fastqe_data.items()):
            for stat_type, emoji_str in sorted(data.items()):
                rows.append(
                    f"<tr>"
                    f'<td style="font-weight:bold;">{s_name}</td>'
                    f"<td>{stat_type}</td>"
                    f'<td style="white-space:nowrap; font-size:1.2em; letter-spacing:2px;">{emoji_str}</td>'
                    f"</tr>"
                )

        if not rows:
            return "<p>No emoji quality data found.</p>"

        return (
            '<table class="table table-striped" style="width:auto;">'
            "<thead><tr>"
            '<th style="width:15%;">Sample</th>'
            '<th style="width:10%;">Statistic</th>'
            '<th style="width:75%;">Quality Emojis</th>'
            "</tr></thead>"
            f"<tbody>{''.join(rows)}</tbody>"
            "</table>"
        )

    def _compute_numeric_stats(self, fastqe_data: dict[str, dict[str, str]]) -> dict[str, dict[str, float]]:
        """Convert emoji quality strings to numeric statistics."""
        numeric_data: dict[str, dict[str, float]] = {}

        for s_name, data in fastqe_data.items():
            # Use 'mean' stat if available, otherwise first available
            emoji_str = data.get("mean") or next(iter(data.values()), None)
            if not emoji_str:
                continue

            scores = self._emoji_to_scores(emoji_str)
            if not scores:
                continue

            numeric_data[s_name] = {
                "avg_quality": sum(scores) / len(scores),
                "min_quality": min(scores),
                "max_quality": max(scores),
                "num_bases": len(scores),
            }

        return numeric_data

    @staticmethod
    def _emoji_to_scores(emoji_str: str) -> list[int]:
        """Convert an emoji quality string to a list of Phred scores."""
        scores = []
        for char in emoji_str:
            if char in EMOJI_TO_PHRED:
                scores.append(EMOJI_TO_PHRED[char])
            # Skip spaces and other non-emoji characters
        return scores
