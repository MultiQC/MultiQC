# https://toolshed.g2.bx.psu.edu/repository?repository_id=13576f42f394cfb6&changeset_revision=c1a33f603da2
# https://github.com/fastqe/fastqe/issues/11
import logging
import os
import re
from typing import Dict, Optional, List

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph
from multiqc.plots.plotly.bar import BarPlotConfig
from multiqc.plots.plotly.line import LinePlotConfig

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    **Note** - MultiQC parses the standard out from FastQE. As such, you must capture the
    FastQE stdout to a file when running to use the MultiQC module.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="FastQE",
            anchor="fastqe",
            href="https://github.com/fastqe/fastqe",
            info="FASTQ sequence quality visualisation with Emoji",
        )
        # Find and load any FastQE reports
        fastqe_data: Dict = dict()
        for f in self.find_log_files("fastqe", filehandles=True):
            self.parse_fastqe_log(f, fastqe_data)

        # Filter to strip out ignored sample names
        fastqe_data = self.ignore_samples(fastqe_data)

        if len(fastqe_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(fastqe_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(fastqe_data, "multiqc_fastqe")

        # Add section with quality scores plot
        self.add_section(
            name="Quality Scores",
            anchor="fastqe-quality",
            description="Quality scores across all bases visualized with emojis",
            plot=self.fastqe_quality_plot(fastqe_data),
        )

    def parse_fastqe_log(self, f, fastqe_data: Dict):
        """Parse the FastQE log output"""
        s_name: Optional[str] = None
        qualities: List[str] = []

        for line in f["f"]:
            # Get sample name from header
            if "Sample Name" in line and "Type" in line and "Qualities" in line:
                continue

            # Parse data line
            fields = line.strip().split("\t")
            if len(fields) == 3:
                s_name = self.clean_s_name(fields[0], f)
                qual_type = fields[1]
                emoji_quals = fields[2]

                if s_name not in fastqe_data:
                    fastqe_data[s_name] = {}
                fastqe_data[s_name][qual_type] = emoji_quals

    @staticmethod
    def fastqe_quality_plot(fastqe_data):
        """Generate quality score visualization"""
        plot_data = []
        for s_name, d in fastqe_data.items():
            if "mean" in d:
                plot_data.append({"name": s_name, "qualities": f'<code style="white-space:nowrap;">{d["mean"]}</code>'})

        return {"id": "fastqe_quality_plot", "plot_type": "html", "data": plot_data}

    # TODO Add to general stats table?
    # TODO Add section?
    # TODO Add test data to test-data repo
