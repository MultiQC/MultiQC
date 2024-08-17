import logging
import re
from typing import Dict, List, Optional

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


def featurecounts_chart(data_keys, data_by_sample):
    """Make the featureCounts assignment rates plot"""

    headers = {}
    for h in data_keys:
        nice_name = h.replace("Unassigned_", "Unassigned: ").replace("_", " ")
        nice_name = re.sub(r"([a-z])([A-Z])", r"\g<1> \g<2>", nice_name)
        headers[h] = {"name": nice_name}

    return bargraph.plot(
        data_by_sample,
        headers,
        {
            "id": "featureCounts_assignment_plot",
            "title": "featureCounts: Assignments",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        },
    )


class MultiqcModule(BaseMultiqcModule):
    """
    As of MultiQC v1.10, the module should also work with output from
    [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html).
    Note that your filenames must end in `.summary` to be discovered.
    See [Module search patterns](#module-search-patterns) for how to customise this.

    Please note that if files are in "Rsubread mode" then lines will be split by any
    whitespace, instead of tab characters. As such, filenames with spaces in will
    cause the parsing to fail.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="featureCounts",
            anchor="featurecounts",
            target="Subread featureCounts",
            href="http://subread.sourceforge.net/",
            info="Counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, "
            "genomic bins and chromosomal locations.",
            doi="10.1093/bioinformatics/btt656",
        )

        data_by_sample: Dict[str, Dict[str, List[int]]] = dict()
        data_keys: List[str] = []
        for f in self.find_log_files("featurecounts"):
            self.parse_featurecounts_report(f, data_by_sample, data_keys)

        # Filter to strip out ignored sample names
        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(data_by_sample, "multiqc_featurecounts")

        # Basic Stats Table
        self.general_stats_table(data_by_sample)

        # Assignment bar plot
        self.add_section(plot=featurecounts_chart(data_keys, data_by_sample))

    def parse_featurecounts_report(self, f, data_by_sample, data_keys):
        """Parse the featureCounts log file."""
        file_names: List[str] = []
        parsed_data: Dict[str, List[int]] = dict()
        split_sep: Optional[str] = "\t"
        for line in f["f"].splitlines():
            this_row: List[int] = list()

            # If this is from Rsubread then the formatting can be very variable
            # Default search pattern is quite generic, so f
            if len(file_names) == 0 and len(line.split(split_sep)) < 2:
                # Split by whitespace and strip quote marks
                # NB: Will break if sample names have whitespace. Rsubread output is so variable that this is difficult to avoid
                split_sep = None

            tokens: List[str] = line.split(split_sep)
            tokens = [token.strip('"') for token in tokens]
            if len(tokens) < 2:
                continue

            key = tokens[0]
            if tokens[0] == "Status":
                for f_name in tokens[1:]:
                    file_names.append(f_name)
            elif len(file_names) + 1 == len(tokens):
                if key not in data_keys:
                    data_keys.append(key)
                for val in tokens[1:]:
                    try:
                        this_row.append(int(val))
                    except ValueError:
                        pass
            if len(this_row) > 0:
                parsed_data[key] = this_row

        # Check that this actually is a featureCounts file, as format and parsing is quite general
        if "Assigned" not in parsed_data.keys():
            return None

        for idx, f_name in enumerate(file_names):
            # Clean up sample name
            s_name = self.clean_s_name(f_name, f)

            # Reorganised parsed data for this sample
            # Collect total count number
            data: Dict[str, float] = dict()
            data["Total"] = 0
            for k in parsed_data:
                data[k] = parsed_data[k][idx]
                data["Total"] += parsed_data[k][idx]

            # Calculate the percent aligned if we can
            try:
                data["percent_assigned"] = (float(data["Assigned"]) / float(data["Total"])) * 100.0
            except (KeyError, ZeroDivisionError):
                pass

            # Add to the main dictionary
            if len(data) > 1:
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name)
                data_by_sample[s_name] = data

    def general_stats_table(self, data_by_sample):
        """Take the parsed stats from the featureCounts report and add them to the
        basic stats table at the top of the report"""

        headers = {
            "Assigned": {
                "title": "Assigned",
                "description": f"Assigned reads ({config.read_count_desc})",
                "scale": "PuBu",
                "shared_key": "read_count",
                "hidden": True,
            },
            "percent_assigned": {
                "title": "Assigned",
                "description": "% Assigned reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
        }
        self.general_stats_addcols(data_by_sample, headers)
