import logging
import os
import re
from typing import Dict, Optional

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.plots.plotly.bar import BarPlotConfig

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    **Note** - MultiQC parses the standard out from Kallisto, _not_ any of its output files
    (`abundance.h5`, `abundance.tsv`, and `run_info.json`). As such, you must capture the
    Kallisto stdout to a file when running to use the MultiQC module.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Kallisto",
            anchor="kallisto",
            href="http://pachterlab.github.io/kallisto/",
            info="Quantifies abundances of transcripts (or more generally, of target sequences) from RNA-Seq data",
            doi="10.1038/nbt.3519",
        )

        # Find and load any Kallisto reports
        kallisto_data: Dict = dict()
        for f in self.find_log_files("kallisto", filehandles=True):
            self.parse_kallisto_log(f, kallisto_data)

        # Filter to strip out ignored sample names
        kallisto_data = self.ignore_samples(kallisto_data)

        if len(kallisto_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(kallisto_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(kallisto_data, "multiqc_kallisto")

        # Basic Stats Table
        self.kallisto_general_stats_table(kallisto_data)

        # Alignment Rate Plot
        self.add_section(plot=self.kallisto_alignment_plot(kallisto_data))

    def parse_kallisto_log(self, f, kallisto_data: Dict):
        s_name: Optional[str] = None
        total_reads: Optional[float] = None
        pseudo_aligned_reads: Optional[float] = None
        frag_length: Optional[float] = None
        for line in f["f"]:
            # Get input filename
            m = re.search(r"\[quant\] will process (pair|file|sample) 1: (\S+)", line)
            if m:
                s_name = self.clean_s_name(os.path.basename(m.group(2)), f)

            if s_name is not None:
                # Alignment rates
                aligned = re.search(r"\[quant\] processed ([\d,]+) reads, ([\d,]+) reads pseudoaligned", line)
                if aligned:
                    total_reads = float(aligned.group(1).replace(",", ""))
                    pseudo_aligned_reads = float(aligned.group(2).replace(",", ""))

                # Paired end fragment lengths
                m = re.search(r"\[quant\] estimated average fragment length: ([\d\.]+)", line)
                if m:
                    frag_length = float(m.group(1).replace(",", ""))

                if "quantifying the abundances" in line:
                    if s_name in kallisto_data:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name)
                    if pseudo_aligned_reads is not None and total_reads is not None:
                        kallisto_data[s_name] = {
                            "total_reads": total_reads,
                            "pseudoaligned_reads": pseudo_aligned_reads,
                            "not_pseudoaligned_reads": total_reads - pseudo_aligned_reads,
                        }
                        if total_reads != 0:
                            kallisto_data[s_name]["percent_aligned"] = (pseudo_aligned_reads / total_reads) * 100
                        else:
                            kallisto_data[s_name]["percent_aligned"] = 0.0
                    if frag_length is not None:
                        kallisto_data[s_name]["fragment_length"] = frag_length
                    s_name = total_reads = pseudo_aligned_reads = frag_length = None

    def kallisto_general_stats_table(self, kallisto_data):
        """Take the parsed stats from the Kallisto report and add it to the
        basic stats table at the top of the report"""

        headers = {
            "fragment_length": {
                "title": "Frag Length",
                "description": "Estimated average fragment length",
                "min": 0,
                "suffix": "bp",
                "scale": "RdYlGn",
            },
            "percent_aligned": {
                "title": "% Aligned",
                "description": "% processed reads that were pseudoaligned",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            },
            "pseudoaligned_reads": {
                "title": f"{config.read_count_prefix} Aligned",
                "description": f"Pseudoaligned reads ({config.read_count_desc})",
                "min": 0,
                "scale": "PuRd",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
        }
        self.general_stats_addcols(kallisto_data, headers)

    @staticmethod
    def kallisto_alignment_plot(kallisto_data):
        # Specify the order of the different possible categories
        keys = {
            "pseudoaligned_reads": {"color": "#437bb1", "name": "Pseudoaligned"},
            "not_pseudoaligned_reads": {"color": "#b1084c", "name": "Not aligned"},
        }

        return bargraph.plot(
            kallisto_data,
            keys,
            BarPlotConfig(
                id="kallisto_alignment",
                title="Kallisto: Alignment Scores",
                ylab="# Reads",
                cpswitch_counts_label="Number of Reads",
            ),
        )
