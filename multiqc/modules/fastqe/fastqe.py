# https://github.com/fastqe/fastqe
# https://github.com/fastqe/fastqe/issues/11
import logging
from typing import Dict, Optional

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
from multiqc.plots.table_object import ColumnDict

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

        # Collect all stat types across samples for column headers
        stat_types = sorted({stat for data in fastqe_data.values() for stat in data})

        headers: Dict[str, ColumnDict] = {}
        for stat in stat_types:
            headers[stat] = {
                "title": stat.capitalize(),
                "description": f"Per-base {stat} quality scores as emoji",
                "scale": False,
                "format": lambda x: f'<span style="font-size: 1.5em; letter-spacing: 2px;">{x}</span>',
            }

        self.add_section(
            name="Quality Scores",
            anchor="fastqe-quality",
            description="Per-base quality scores represented as emoji. Each emoji maps to a Phred quality score range.",
            helptext="""
Each emoji maps to a Phred quality score. The default scale from 0 to 41 is:

🚫 ❌ 👺 💔 🙅 👾 👿 💀 👻 🙈 🙉 🙊 🐵 😿 😾 🙀 💣 🔥 😡 💩 🚨 😀 😅 😏 😊 😙 😗 😚 😃 😘 😆 😄 😋 😝 😛 😜 😉 😁 😎 😍 🤩 😍

The binned scale (`--bin`) compresses scores into 8 ranges:
🚫 (0-1), 💀 (2-9), 💩 (10-19), ⚠️ (20-24), 😄 (25-29), 😆 (30-34), 😎 (35-39), 😍 (40+)

FastQE supports `--mean` (default), `--min`, and `--max` statistics.
Each will appear as a separate column in the table.
            """,
            plot=table.plot(
                fastqe_data,
                headers,
                pconfig={
                    "namespace": "FastQE",
                    "id": "fastqe_quality_table",
                    "title": "FastQE: Quality Scores",
                    "parse_numeric": False,
                    "sort_rows": False,
                    "no_violin": True,
                },
            ),
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
