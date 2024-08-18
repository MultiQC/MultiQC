import logging
import os
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Note that the versions < 1.7.8 use the basename of the file path to distinguish samples, whereas newer
    versions produce logfiles with a sample identifer that gets parsed by MultiQC.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ClipAndMerge",
            anchor="clipandmerge",
            href="http://www.github.com/apeltzer/ClipAndMerge",
            info="Adapter clipping and read merging for ancient DNA data.",
            doi="10.1186/s13059-016-0918-z",
        )

        # Find and load any ClipAndMerge reports
        self.clipandmerge_data = dict()
        for f in self.find_log_files("clipandmerge"):
            self.parse_clipandmerge_log(f)

        # Filter to strip out ignored sample names
        self.clipandmerge_data = self.ignore_samples(self.clipandmerge_data)

        if len(self.clipandmerge_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.clipandmerge_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.clipandmerge_data, "multiqc_clipandmerge")

        # Basic Stats Table
        self.clipandmerge_general_stats_table()

        # Alignment Rate Plot
        self.add_section(plot=self.clipandmerge_alignment_plot())

    def parse_clipandmerge_log(self, f):
        regexes = {
            "usable_reads": r"Number of usable reads in the output file\(s\):\s+(\d+)",
            "merged_reads": r"Number of usable merged reads:\s+(\d+)",
            "percentage": r"Percentage of usable merged reads:\s+(\d+\.\d+)",
            "usable_not_merged_forward": r"Number of usable not merged forward reads:\s+(\d+)",
            "usable_not_merged_reverse": r"Number of usable not merged reverse reads:\s+(\d+)",
            "usable_forward_no_pairing_reverse": r"Number of usable forward reads with no pairing reverse read:\s+(\d+)",
            "usable_reverse_no_pairing_forward": r"Number of usable reverse reads with no pairing forward read:\s+(\d+)",
            "identifier": r"SampleID:\s+(\S+)",
            "version": r"ClipAndMerge \(v. ([\d\.]+)\)",
        }

        parsed_data = dict()
        for k, r in regexes.items():
            r_search = re.search(r, f["f"], re.MULTILINE)
            if r_search:
                if k == "versions":
                    parsed_data[k] = r_search.group(1)
                    continue

                try:
                    parsed_data[k] = float(r_search.group(1))
                except ValueError:
                    parsed_data[k] = r_search.group(1)

        if len(parsed_data) > 0:
            s_name = self.clean_s_name(os.path.basename(f["root"]), f)
            if "identifier" in parsed_data:
                s_name = self.clean_s_name(parsed_data["identifier"], f)
            self.clipandmerge_data[s_name] = parsed_data

            if "version" in parsed_data:
                self.add_software_version(parsed_data["version"], s_name)

            self.add_data_source(f)

    def clipandmerge_general_stats_table(self):
        """Take the parsed stats from the ClipAndMerge report and add it to the
        basic stats table at the top of the report"""

        headers = {
            "percentage": {
                "title": "% Merged",
                "description": "Percentage of reads merged",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Greens",
                "format": "{:,.2f}",
            }
        }
        self.general_stats_addcols(self.clipandmerge_data, headers)

    def clipandmerge_alignment_plot(self):
        # Specify the order of the different possible categories
        keys = {
            "merged_reads": {"name": "Merged Reads"},
            "usable_not_merged_forward": {"name": "Usable, not merged (forward)"},
            "usable_not_merged_reverse": {"name": "Usable, not merged (reverse)"},
            "usable_forward_no_pairing_reverse": {"name": "Usable forward-only"},
            "usable_reverse_no_pairing_forward": {"name": "Usable reverse-only"},
        }

        # Config for the plot
        config = {
            "id": "clipandmerge_rates",
            "title": "ClipAndMerge: Read merging results",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
            "hide_empty": False,
        }

        return bargraph.plot(self.clipandmerge_data, keys, config)
