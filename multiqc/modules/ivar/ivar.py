import logging
import re

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import heatmap

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This module parses the output from the `ivar trim` command and creates a table view.
    Both output from V1 and V2 of the tool are supported and parsed accordingly.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="iVar",
            anchor="iVar",
            href="https://github.com/andersen-lab/ivar",
            info="Functions for viral amplicon-based sequencing.",
            doi="10.1101/383513",
        )

        # Find and load iVar trim results
        self.ivar_data = dict()
        self.ivar_primers = dict()
        for f in self.find_log_files("ivar/trim", filehandles=True):
            self.parse_ivar(f)

        # Filter to strip out ignored sample names
        self.ivar_data = self.ignore_samples(self.ivar_data)
        self.ivar_primers = self.ignore_samples(self.ivar_primers)

        # Stop if no files are found
        num_samples = max(len(self.ivar_data), len(self.ivar_primers))
        if num_samples == 0:
            raise ModuleNoSamplesFound

        # Log number of reports
        log.info(f"Found {num_samples} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed data to a file
        self.write_data_file(self.ivar_data, "multiqc_ivar_summary")
        self.write_data_file(self.ivar_primers, "multiqc_ivar_primers")

        # Basic Stats Table
        self.ivar_general_stats_table()

        # Heatmap info
        self.primer_heatmap()

    # Parse an ivar report
    def parse_ivar(self, f):
        """Parse iVar log file, either v1.1_beta or v1.2"""
        count_regexes = {
            "mapped_reads": re.compile(r"(?:Found\s)(\d+)(?:\smapped)"),
            "trimmed_reads": re.compile(r"(?:Trimmed primers from )(?:\d+\.\d+\% \()?(\d+)"),
            "reads_outside_primer_region": re.compile(
                r"^(?:\d+\.\d+\% \()?(\d+)(?:\))?(?:.*[of]?)reads\s(?:that\s)?started"
            ),
            "reads_too_short_after_trimming": re.compile(
                r"^(?:\d+\.\d+\% \()?(\d+)(?:\))?(?:.*[of]?)reads\swere(?: quality trimmed | shortened)"
            ),
        }
        primer_regex = re.compile(r"^(.*)(?:\t+)(\d+$)")
        parsed_data = dict()
        primers = dict()
        for line in f["f"]:
            # Search count regexes for stats
            for k, count_regex in count_regexes.items():
                count_match = count_regex.search(line)
                if count_match:
                    parsed_data[k] = int(count_match.group(1))

            # Try to match the primer regex
            primer_match = primer_regex.search(line)
            if primer_match:
                primers[primer_match.group(1)] = int(primer_match.group(2))

        if parsed_data is not None and len(parsed_data) > 0:
            if f["s_name"] in self.ivar_data:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.ivar_data[f["s_name"]] = parsed_data
            self.add_data_source(f, section="trimming")

        if primers is not None and len(primers) > 0:
            self.ivar_primers[f["s_name"]] = primers
            self.add_data_source(f, section="trimming")

    # Add to general stats table
    def ivar_general_stats_table(self):
        """Take the parsed stats from the iVAR report and add it to the
        basic stats table"""

        headers = {
            "reads_too_short_after_trimming": {
                "title": f"{config.read_count_prefix} Too short",
                "description": "Number of reads too short (<30bp) after primer trimming ({})".format(
                    config.read_count_desc
                ),
                "scale": "OrRd",
                "shared_key": "read_counts",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "reads_outside_primer_region": {
                "title": f"{config.read_count_prefix} Outside primer",
                "description": f"Number of reads outside the primer region ({config.read_count_desc})",
                "scale": "YlOrBr",
                "shared_key": "read_counts",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "trimmed_reads": {
                "title": f"{config.read_count_prefix} Primer trimmed",
                "description": "Total number of reads where primer trimming was performed. ({})".format(
                    config.read_count_desc
                ),
                "scale": "Purples",
                "shared_key": "read_counts",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "mapped_reads": {
                "title": f"{config.read_count_prefix} Mapped",
                "description": f"Total number of mapped reads in iVar input. ({config.read_count_desc})",
                "scale": "PuBu",
                "shared_key": "read_counts",
                "modify": lambda x: x * config.read_count_multiplier,
            },
        }
        self.general_stats_addcols(self.ivar_data, headers)

    def primer_heatmap(self):
        """Heatmap showing information on each primer found for every sample"""
        # Top level dict contains sample IDs + dict(primer, counts)

        final_data = []
        final_xcats = []
        final_ycats = []

        for i, (k, v) in enumerate(self.ivar_primers.items()):
            final_ycats.append(k)
            tmp_prim_val = []
            for prim, val in v.items():
                if i == 0:
                    final_xcats.append(prim)
                tmp_prim_val.append(val)
            final_data.append(tmp_prim_val)

        if len(final_xcats) != len(final_data[0]):
            log.error(f"Number of primers do not match: {len(final_xcats)} vs {len(final_data[0])}")

        if len(final_ycats) != len(final_data):
            log.error(f"Number of samples and primers do not match: {len(final_ycats)} vs {len(final_data)}")

        if self.ivar_primers is not None:
            pconfig = {
                "id": "ivar-primer-count-heatmap",
                "title": "iVar: Number of primers found for each sample",
                "tt_decimals": 0,
                "square": False,
                "xcats_samples": False,
            }

            self.add_section(
                name="iVar Primer Counts",
                anchor="ivar-primers-heatmap",
                description="Counts observed for each primer per sample.",
                helptext="This lists the number of times a specific primer was found in the respective sample.",
                plot=heatmap.plot(final_data, final_xcats, final_ycats, pconfig),
            )
