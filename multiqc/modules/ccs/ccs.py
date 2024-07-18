import json
import logging
import re
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="CCS",
            anchor="ccs",
            href="https://github.com/PacificBiosciences/ccs",
            info="PacBio tool that generates highly accurate single-molecule consensus reads (HiFi Reads)",
            extra="""
            CCS takes multiple subreads of the same SMRTbell molecule and combines them
            using a statistical model to produce one highly accurate consensus sequence,
            also called HiFi read, with base quality values. This tool powers the Circular
            Consensus Sequencing workflow in SMRT Link.
            """,
            # Can't find a DOI // doi=
        )

        # To store the mod data
        self.ccs_data: Dict = dict()
        self.parse_v4_log_files()
        self.parse_v5_log_files()
        self.ccs_data = self.ignore_samples(self.ccs_data)

        # If we found no data
        if not self.ccs_data:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.ccs_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.write_data_file(self.ccs_data, "multiqc_ccs_report")
        self.add_general_stats()
        self.add_sections()

    def parse_v4_log_files(self):
        for f in self.find_log_files("ccs/v4", filehandles=True):
            data = parse_PacBio_log(f["f"])
            v5_data = self.convert_to_v5(data)
            self.ccs_data[f["s_name"]] = v5_data
            self.add_data_source(f)

    def parse_v5_log_files(self):
        for f in self.find_log_files("ccs/v5", filehandles=True):
            v5_data = json.load(f["f"])
            self.ccs_data[f["s_name"]] = v5_data
            self.add_data_source(f)

    def add_general_stats(self):
        gstats_data = {}
        for s_name, attrs in self.ccs_data.items():
            gstats_data[s_name] = {}
            for attr in attrs["attributes"]:
                if "zmw_input" in attr["id"]:
                    gstats_data[s_name]["zmw_input"] = attr["value"]
                if "zmw_passed_yield" in attr["id"]:
                    gstats_data[s_name]["zmw_passed_yield"] = attr["value"]
            try:
                gstats_data[s_name]["zmw_pct_passed_yield"] = (
                    gstats_data[s_name]["zmw_passed_yield"] / gstats_data[s_name]["zmw_input"]
                ) * 100.0
            except (KeyError, ZeroDivisionError):
                pass

        headers = {
            "zmw_pct_passed_yield": {
                "title": "ZMWs %PF",
                "description": "Percent of ZMWs passing filters",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn",
                "suffix": "%",
            },
            "zmw_passed_yield": {
                "title": "ZMWs PF",
                "description": "ZMWs pass filters",
                "scale": "BuGn",
                "format": "{:,d}",
                "shared_key": "zmw_count",
                "hidden": True,
            },
            "zmw_input": {
                "title": "ZMWs input",
                "description": "ZMWs input",
                "scale": "Purples",
                "format": "{:,d}",
                "shared_key": "zmw_count",
            },
        }
        self.general_stats_addcols(gstats_data, headers)

    def add_sections(self):
        # First we gather all the filters we encountered
        all_filters = dict()
        for s_name in self.ccs_data:
            for filter_reason in self.filter_and_pass(self.ccs_data[s_name]):
                all_filters[filter_reason] = 0

        # Then we add the counts for each filter to the plot data
        plot_data = dict()
        for s_name, data in self.ccs_data.items():
            plot_data[s_name] = dict()
            for reason in all_filters:
                for attribute in data["attributes"]:
                    if attribute["name"] == reason:
                        plot_data[s_name][reason] = attribute["value"]
                        all_filters[reason] += attribute["value"]
                        break
                # If we didn't find it, set to zero
                else:
                    plot_data[s_name][reason] = 0

        # Sort the categories based on the total counts
        # Smallest first, for Log10 view
        plot_cats = []
        for reason, count in sorted(all_filters.items(), key=lambda item: item[1], reverse=True):
            plot_cats.append(reason)

        # Plot configuration
        config = {
            "id": "ccs-filter-graph",
            "title": "CCS: ZMW results",
            "ylab": "Number of ZMWs",
            "xlab": "CCS report file",
            "logswitch": True,
        }

        self.add_section(
            name="ZMWs filtered and passed",
            anchor="ccs-filter",
            description="The number of ZMWs that failed or passed all filters for each of the subreads files.",
            helptext="""
                The number of ZMWs that passed all filters is shown as **ZMWs pass filters**. All other categories that
                are shown in the graph represent the number of ZMWs that were dropped for the specified reason.
            """,
            plot=bargraph.plot(plot_data, plot_cats, config),
        )

    def filter_and_pass(self, data):
        """Gather the reasons why ZMWs were failed or passed"""
        reasons = dict()

        # We only have to use the attributes
        attributes = data["attributes"]

        # Add filtere reasons (id starts with filtered_) and total passed
        for entry in attributes:
            if (
                entry["id"].startswith("filtered")
                or entry["id"] == "zmw_passed_yield"
                or entry["id"].startswith("ccs_processing.filtered")
                or entry["id"] == "ccs_processing.zmw_passed_yield"
            ):
                reasons[entry["name"]] = entry["value"]

        return reasons

    def convert_to_v5(self, data):
        """Convert the v4 format to the new CCS v5 json format"""
        # Initialise the v5 format dictionary
        v5 = dict()
        v5["id"] = "ccs_processing"
        attributes = list()
        v5["attributes"] = attributes

        # Update names for top level entries, they have been changed in v5
        data["ZMWs pass filters"] = data["ZMWs generating CCS"]
        data["ZMWs fail filters"] = data["ZMWs filtered"]

        # Add the top level values to the attributes
        for name in ("ZMWs input", "ZMWs pass filters", "ZMWs fail filters"):
            count = data[name]["count"]
            attributes.append(self.make_v5_entry(name, count))

        # Add the reasons for filtering to the attributes
        for reason in data["ZMWs filtered"]["Exclusive ZMW counts"]:
            count = data["ZMWs filtered"]["Exclusive ZMW counts"][reason]["count"]
            attributes.append(self.make_v5_entry(reason, count))

        return v5

    def make_v5_entry(self, name, count):
        """Make a v5 output entry based on name and count"""
        # Dictionary to map the report names to the v5 id annotation
        name_to_id = {
            "ZMWs input": "zmw_input",
            "ZMWs pass filters": "zmw_passed_yield",
            "ZMWs fail filters": "zmw_filtered_yield",
            "Below SNR threshold": "filtered_poor_snr",
            "Median length filter": "filtered_median_length_filter",
            "Lacking full passes": "filtered_too_few_passes",
            "Heteroduplex insertions": "filtered_heteroduplexes",
            "Coverage drops": "filtered_coverage_drops",
            "Insufficient draft cov": "filtered_insufficient_draft_cov",
            "Draft too different": "filtered_draft_too_different",
            "Draft generation error": "filtered_draft_failure",
            "Draft above --max-length": "filtered_draft_too_long",
            "Draft below --min-length": "filtered_draft_too_short",
            "Reads failed polishing": "filtered_read_failed_polish",
            "Empty coverage windows": "filtered_empty_window_during_polishing",
            "CCS did not converge": "filtered_non_convergent",
            "CCS below minimum RQ": "filtered_below_rq",
            "Unknown error": "filtered_unknown_error",
        }

        return {"name": name, "id": name_to_id[name], "value": count}


def parse_PacBio_log(file_content):
    """Parse ccs log file"""
    data = dict()
    # This is a local dictionary to store which annotations belong to which
    # result dictionary. This will be used to structure the output, but will
    # not be part of the output itself
    annotations = dict()
    current_annotation = None

    for line in file_content:
        # Get rid of trailing newlines
        line = line.strip("\n")
        # Did we enter a new section with annotations for an earlier result?
        # If so, we will only add an empty dictionary with the correct name
        # These field are of the format "something something for (A):"
        section_header_pattern = " for [(][A-Z][)][:]$"
        if re.search(section_header_pattern, line):
            linedata = parse_line(line)
            ann = linedata["annotation"]
            # Cut off the ' for (B):' part
            name = line[:-9]
            # We make a new heading with the current name under the data that
            # matches the current annotation
            current_annotation = dict()
            # We keep the dictonary accessible under 'current_annotation',
            # so we can keep adding new data to it without having to keep track
            # of where it belongs
            annotations[ann][name] = current_annotation
            continue

        linedata = parse_line(line)

        # If we got no data, we reached the end of the annotated section
        if not linedata:
            current_annotation = None
            continue

        # Lets get the name of the data
        name = linedata.pop("name")
        # If we are in an annotated section, we add the data to the current
        # annotation
        if current_annotation is not None:
            current_annotation[name] = linedata
        # Otherwise, we add the newfound annotation to the dictionary in case
        # we find a corresponding section later on.
        # The data from the current line we add directly to the output data
        else:
            annotation = linedata.pop("annotation")
            annotations[annotation] = linedata
            data[name] = linedata

    return data


def parse_line(line):
    """Parse a line from the ccs log file"""
    data = dict()

    # If we got an empty line to parse
    if not line.strip():
        return data

    # Split the line on the colon character
    key, values = line.strip().split(":")

    # The key can have multiple parts
    keys = key.strip().split()

    # Does the key have an annotation (A), (B) etc
    if re.fullmatch("[(][A-Z][)]", keys[-1]):
        # We store the annotation without the bracets
        data["annotation"] = keys[-1][1:-1]
        # And we add the rest of the key as name
        data["name"] = " ".join(keys[:-1])
    # Otherwise, we just store the name
    else:
        data["name"] = " ".join(keys)

    # Parsing the values
    values = values.strip().split()
    # Are there are no values we are done
    if not values:
        return data

    # If there is a single value
    if len(values) == 1:
        value = values.pop()
        # Is the value a percentage
        if value.endswith("%"):
            data["percentage"] = float(value[:-1])
        # Otherwise, it is an integer
        else:
            data["count"] = int(value)
    elif len(values) == 2:
        # If there are multiple values, they are in the format
        # count (percentage%)
        count = values[0]
        percentage = values[1]

        # The percentage should be in the format:
        # (12.34%) or (100%) or (0%) or (-nan%)
        # So we can remove the bracets first
        percentage = percentage[1:-1]
        # Then we make sure it is one of the supported patterns
        assert re.fullmatch(r"(\d+\.\d+%|\d+%|-nan%)", percentage)

        # If the percentage is this weird nan, set it to 0
        if percentage == "-nan%":
            data["percentage"] = 0.0
        else:  # Otherwise, we cut of the % sign and convert to float
            data["percentage"] = float(percentage[:-1])
        # Add the count as an integer
        data["count"] = int(count)

    return data
