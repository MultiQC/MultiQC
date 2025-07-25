import logging
import re
from typing import Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.types import Anchor, ColumnKey

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses the stderr output, that can be captured by directing it to a file e.g.:

    ```sh
    trimmomatic command 2> trim_out.log
    ```

    By default, the module generates the sample names based on the input FastQ file names in
    the command line used by Trimmomatic. If you prefer, you can tell the module to use
    the filenames as sample names instead. To do so, use the following config option:

    ```yaml
    use_filename_as_sample_name: true
    ```

    Note: The old `trimmomatic.s_name_filenames` option is deprecated and will be removed in a future version.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Trimmomatic",
            anchor=Anchor("trimmomatic"),
            href="http://www.usadellab.org/cms/?page=trimmomatic",
            info="Read trimming tool for Illumina NGS data.",
            doi="10.1093/bioinformatics/btu170",
        )

        # Parse logs
        self.trimmomatic: Dict = dict()
        for f in self.find_log_files("trimmomatic", filehandles=True):
            self.parse_trimmomatic(f)
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.trimmomatic = self.ignore_samples(self.trimmomatic)
        if len(self.trimmomatic) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.trimmomatic)} logs")

        self.write_data_file(self.trimmomatic, "multiqc_trimmomatic")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add drop rate to the general stats table
        self.general_stats_addcols(
            self.trimmomatic,
            {
                ColumnKey("dropped_pct"): {
                    "title": "% Dropped",
                    "description": "% Dropped reads",
                    "max": 100,
                    "min": 0,
                    "suffix": "%",
                    "scale": "OrRd",
                }
            },
        )

        # Make barplot
        self.trimmomatic_barplot()

    def parse_trimmomatic(self, f):
        s_name = None
        should_use_filename = False
        if isinstance(config.use_filename_as_sample_name, list):
            # Check for module anchor
            if self.anchor in config.use_filename_as_sample_name:
                should_use_filename = True
        elif config.use_filename_as_sample_name is True:
            should_use_filename = True
        elif getattr(config, "trimmomatic", {}).get("s_name_filenames", False):
            # Deprecated option - warn user
            log.warning(
                "The 'trimmomatic.s_name_filenames' config option is deprecated. Use the global 'use_filename_as_sample_name' option instead."
            )
            should_use_filename = True

        if should_use_filename:
            s_name = f["s_name"]
        for line in f["f"]:
            # Get the sample name
            if s_name is None and line.startswith(
                tuple(f"Trimmomatic{x}E: Started with arguments:" for x in ["S", "P"])
            ):
                is_pe = line.startswith("TrimmomaticPE")
                args = line.strip().split()
                FQ_EXTS = ".fastq", ".fq", ".gz", ".dat"
                if not any(x.endswith(FQ_EXTS) for x in args):
                    # Try looking on the next line instead, sometimes have a line break (see issue #212)
                    line = next(f["f"])
                    args = line.strip().split()
                if any(x.endswith(FQ_EXTS) for x in args):
                    # For PE, first two fastq files are the input paths; for SE, it's just the first one
                    input_paths = [x for x in args if x.endswith(FQ_EXTS)][: 2 if is_pe else 1]
                    s_name = self.clean_s_name(input_paths, f)
                    if s_name in self.trimmomatic:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")

            # Get single end stats
            if "Input Reads" in line and s_name is not None:
                match = re.search(
                    r"Input Reads: (\d+) Surviving: (\d+) \(([\d\.,]+)%\) Dropped: (\d+) \(([\d\.,]+)%\)", line
                )
                if match:
                    self.trimmomatic[s_name] = {
                        "input_reads": float(match.group(1)),
                        "surviving": float(match.group(2)),
                        "surviving_pct": float(match.group(3).replace(",", ".")),
                        "dropped": float(match.group(4)),
                        "dropped_pct": float(match.group(5).replace(",", ".")),
                    }
                    s_name = None

            # Get paired end stats
            if "Input Read Pairs" in line and s_name is not None:
                match = re.search(
                    r"Input Read Pairs: (\d+) Both Surviving: (\d+) \(([\d\.,]+)%\) Forward Only Surviving: (\d+) \(([\d\.,]+)%\) Reverse Only Surviving: (\d+) \(([\d\.,]+)%\) Dropped: (\d+) \(([\d\.,]+)%\)",
                    line,
                )
                if match:
                    self.trimmomatic[s_name] = {
                        "input_read_pairs": float(match.group(1)),
                        "surviving": float(match.group(2)),
                        "surviving_pct": float(match.group(3).replace(",", ".")),
                        "forward_only_surviving": float(match.group(4)),
                        "forward_only_surviving_pct": float(match.group(5).replace(",", ".")),
                        "reverse_only_surviving": float(match.group(6)),
                        "reverse_only_surviving_pct": float(match.group(7).replace(",", ".")),
                        "dropped": float(match.group(8)),
                        "dropped_pct": float(match.group(9).replace(",", ".")),
                    }
                    s_name = None

    def trimmomatic_barplot(self):
        # Specify the order of the different possible categories
        keys = {
            "surviving": {"color": "#437bb1", "name": "Surviving Reads"},
            "both_surviving": {"color": "#f7a35c", "name": "Both Surviving"},
            "forward_only_surviving": {"color": "#e63491", "name": "Forward Only Surviving"},
            "reverse_only_surviving": {"color": "#b1084c", "name": "Reverse Only Surviving"},
            "dropped": {"color": "#7f0000", "name": "Dropped"},
        }

        # Config for the plot
        pconfig = {
            "id": "trimmomatic_plot",
            "title": "Trimmomatic: Surviving Reads",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        self.add_section(plot=bargraph.plot(self.trimmomatic, keys, pconfig))
