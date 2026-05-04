# https://github.com/fastqe/fastqe
# https://github.com/fastqe/fastqe/issues/11
import html
import logging
from typing import Dict, Optional

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    FastQE uses emoji to represent FASTQ sequence quality scores, providing a
    fun and visually intuitive way to assess sequencing data quality.

    The module parses the tab-separated output from FastQE and displays the
    emoji quality strings for each sample in the report.

    **Note** — MultiQC parses the standard output from FastQE. You must capture
    FastQE stdout to a file when running, for example:

    ```bash
    fastqe input.fastq > fastqe_output.txt
    ```

    The saved file must have `fastqe` somewhere in the file name.
    """

    def __init__(self):
        super().__init__(
            name="FastQE",
            anchor="fastqe",
            href="https://github.com/fastqe/fastqe",
            info="Uses emoji to represent FASTQ sequence quality scores.",
        )

        fastqe_data: Dict[str, Dict[str, str]] = {}
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
            description="Per-base quality scores represented as emoji. Each emoji maps to a Phred quality score range.",
            plot=self._fastqe_emoji_table(fastqe_data),
        )

        self.write_data_file(fastqe_data, "multiqc_fastqe")

    def _parse_fastqe_log(self, f) -> Optional[Dict[str, Dict[str, str]]]:
        """Parse FastQE TSV: Filename\\tStatistic\\tQualities -> {sample: {stat: emoji}}"""
        data: Dict[str, Dict[str, str]] = {}

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
    def _fastqe_emoji_table(fastqe_data: Dict[str, Dict[str, str]]) -> str:
        """Generate an HTML table showing emoji quality strings per sample."""
        rows = []
        for s_name, data in sorted(fastqe_data.items()):
            for stat_type, emoji_str in sorted(data.items()):
                rows.append(
                    f"<tr>"
                    f'<td style="font-weight:bold;">{html.escape(s_name)}</td>'
                    f"<td>{html.escape(stat_type)}</td>"
                    f'<td style="white-space:nowrap; font-size:1.2em; letter-spacing:2px;">{html.escape(emoji_str)}</td>'
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
