import json
import logging
import re
from html import escape

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)

# Attribute ID prefixes for the ZMW counts (default) and for the double- and
# single-strand read counts reported with `ccs --hd-finder`
GROUPS = {"zmw": "ZMWs", "ds_reads": "Double-strand reads", "ss_reads": "Single-strand reads"}


class MultiqcModule(BaseMultiqcModule):
    """
    Supports the `ccs` report from CCS v4 (text) and from v5 onwards (JSON, `--report-json`),
    including v8. Reports from `ccs --hd-finder` count double-strand and single-strand reads
    instead of ZMWs, and these are plotted separately.

    CCS v8 reports also include HiFi (Q20+), <Q20 and Q30+ read metrics, which are shown in a
    table and a stacked bar plot when present.
    """

    def __init__(self):
        super().__init__(
            name="CCS",
            anchor="ccs",
            href="https://ccs.how/",
            info="PacBio tool that generates highly accurate single-molecule consensus reads (HiFi Reads)",
            extra="""
            CCS takes multiple subreads of the same SMRTbell molecule and combines them
            using a statistical model to produce one highly accurate consensus sequence,
            also called HiFi read, with base quality values. This tool powers the Circular
            Consensus Sequencing workflow in SMRT Link.
            """,
            # Can't find a DOI // doi=
            license="BSD 3-Clause Clear License",
            license_url="https://github.com/PacificBiosciences/ccs/blob/master/LICENSE",
        )

        # To store the mod data
        self.ccs_data: dict = {}
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

        self.parsed = {s_name: split_attributes(d["attributes"]) for s_name, d in self.ccs_data.items()}
        self.add_general_stats()
        self.add_filter_section()
        self.add_hifi_sections()

        self.write_data_file(self.ccs_data, "multiqc_ccs_report")

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
        for s_name, (groups, other) in self.parsed.items():
            gstats_data[s_name] = {}
            if "zmw" in groups:
                gstats_data[s_name]["zmw_input"] = groups["zmw"]["input"]["value"]
                gstats_data[s_name]["zmw_passed_yield"] = groups["zmw"]["passed_yield"]["value"]
            # --hd-finder: the two strand groups together make up all reads
            strand_groups = [groups[g] for g in ("ds_reads", "ss_reads") if g in groups]
            if strand_groups:
                gstats_data[s_name]["reads_input"] = sum(g["input"]["value"] for g in strand_groups)
                gstats_data[s_name]["reads_passed_yield"] = sum(g["passed_yield"]["value"] for g in strand_groups)
            for unit in ("zmw", "reads"):
                if gstats_data[s_name].get(f"{unit}_input"):
                    gstats_data[s_name][f"{unit}_pct_passed_yield"] = (
                        gstats_data[s_name][f"{unit}_passed_yield"] / gstats_data[s_name][f"{unit}_input"] * 100.0
                    )
            for key in ("number_of_ccs_reads_q20", "total_number_of_ccs_bases_q20", "median_qv_q20"):
                if key in other:
                    gstats_data[s_name][key] = other[key]["value"]

        headers = {}
        for unit, label in (("zmw", "ZMWs"), ("reads", "Reads")):
            headers[f"{unit}_pct_passed_yield"] = {
                "title": f"{label} %PF",
                "description": f"Percent of {label.lower()} passing filters",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn",
                "suffix": "%",
            }
            headers[f"{unit}_passed_yield"] = {
                "title": f"{label} PF",
                "description": f"{label} pass filters",
                "scale": "BuGn",
                "format": "{:,d}",
                "shared_key": f"{unit}_count",
                "hidden": True,
            }
            headers[f"{unit}_input"] = {
                "title": f"{label} input",
                "description": f"{label} input",
                "scale": "Purples",
                "format": "{:,d}",
                "shared_key": f"{unit}_count",
            }
        headers["reads_pct_passed_yield"]["description"] += " (double- and single-strand reads, --hd-finder)"
        headers["number_of_ccs_reads_q20"] = {
            "title": "HiFi reads",
            "description": "Number of HiFi (Q20+) reads",
            "scale": "Blues",
            "format": "{:,d}",
        }
        headers["total_number_of_ccs_bases_q20"] = {
            "title": "HiFi yield",
            "description": "HiFi (Q20+) yield (bp)",
            "scale": "Greens",
            "format": "{:,d}",
            "hidden": True,
        }
        headers["median_qv_q20"] = {
            "title": "HiFi QV",
            "description": "Median read quality of HiFi (Q20+) reads",
            "min": 0,
            "max": 40,
            "scale": "RdYlGn",
            "hidden": True,
        }
        self.general_stats_addcols(gstats_data, headers)

    def add_filter_section(self):
        plot_data: dict[str, dict] = {group: {} for group in GROUPS}
        cats: dict[str, dict] = {}
        totals: dict[str, int] = {}
        for s_name, (groups, _) in self.parsed.items():
            for group, attrs in groups.items():
                for key, attr in attrs.items():
                    if key != "passed_yield" and (not key.startswith("filtered_") or key == "filtered_yield"):
                        continue
                    plot_data[group].setdefault(s_name, {})[key] = attr["value"]
                    # --hd-finder names are prefixed with "DS reads: " / "SS reads: "
                    cats.setdefault(key, {"name": attr["name"].split(": ")[-1]})
                    totals[key] = totals.get(key, 0) + attr["value"]
        # Passed first, then filter reasons largest first
        cats = {key: cats[key] for key in sorted(cats, key=lambda k: (k != "passed_yield", -totals[k]))}
        if "passed_yield" in cats:
            cats["passed_yield"].update({"name": "Pass filters", "color": "#2ca02c"})
        present = [group for group in GROUPS if plot_data[group]]

        self.add_section(
            name="ZMWs filtered and passed",
            anchor="ccs-filter",
            description="The number of ZMWs (or reads, with `--hd-finder`) that failed or passed all filters for each report.",
            helptext="""
                The number that passed all filters is shown as **Pass filters**. All other categories that
                are shown in the graph represent the number of ZMWs (or reads, with `--hd-finder`) that were dropped for the specified reason.

                When CCS is run with `--hd-finder`, it reports double-strand and single-strand reads instead
                of ZMWs. These are shown as separate datasets.
            """,
            plot=bargraph.plot(
                [plot_data[group] for group in present],
                cats,
                {
                    "id": "ccs-filter-graph",
                    "title": "CCS: Filter results",
                    "ylab": "Number of ZMWs",
                    "xlab": "CCS report file",
                    "logswitch": True,
                    "data_labels": [
                        {"name": GROUPS[g], "ylab": "Number of ZMWs" if g == "zmw" else "Number of reads"}
                        for g in present
                    ],
                },
            ),
        )

    def add_hifi_sections(self):
        """Read metrics by quality bin, reported from CCS v8"""
        metrics = {s_name: other for s_name, (_, other) in self.parsed.items() if "number_of_ccs_reads_q20" in other}
        if not metrics:
            return

        visible = {
            "number_of_ccs_reads_q20",
            "total_number_of_ccs_bases_q20",
            "mean_ccs_readlength_q20",
            "median_qv_q20",
            "mean_npasses_q20",
            "base_percentage_q30",
        }
        scales = ["Blues", "Greens", "Purples", "Oranges", "BuPu", "YlGn"]
        headers: dict[str, dict] = {}
        for other in metrics.values():
            for key, attr in other.items():
                # Titles are rendered as HTML but also used as column names in the exported data file
                title = attr["name"].replace("<Q20", "Below Q20").replace(">=Q30", "Q30+")
                headers.setdefault(
                    key,
                    {
                        "title": escape(title),
                        "format": "{:,.0f}",
                        "scale": scales[len(headers) % len(scales)],
                        "hidden": key not in visible,
                    },
                )
                if attr["name"].endswith("(%)"):
                    headers[key].update({"suffix": "%", "min": 0, "max": 100, "format": "{:,.1f}"})
        self.add_section(
            name="Read metrics",
            anchor="ccs-read-metrics",
            description="Read count, yield, length and quality for HiFi (Q20+), <Q20 and Q30+ reads.",
            plot=table.plot(
                {s_name: {k: a["value"] for k, a in other.items()} for s_name, other in metrics.items()},
                headers,
                {"id": "ccs-read-metrics-table", "title": "CCS: Read metrics"},
            ),
        )

        # Q30+ reads are a subset of HiFi (Q20+) reads
        bins_data = []
        for unit in ("number_of_ccs_reads", "total_number_of_ccs_bases"):
            bins_data.append(
                {
                    s_name: {
                        "q30": other[f"{unit}_q30"]["value"],
                        "q20": other[f"{unit}_q20"]["value"] - other[f"{unit}_q30"]["value"],
                        "lq": other[f"{unit}_lq"]["value"],
                    }
                    for s_name, other in metrics.items()
                }
            )
        self.add_section(
            name="Read quality",
            anchor="ccs-read-quality",
            description="Reads and bases in each read quality bin. HiFi reads are those of Q20 and above.",
            plot=bargraph.plot(
                bins_data,
                {
                    "q30": {"name": "Q30+", "color": "#2ca02c"},
                    "q20": {"name": "Q20 to Q30", "color": "#98df8a"},
                    "lq": {"name": "<Q20", "color": "#d62728"},
                },
                {
                    "id": "ccs-read-quality-graph",
                    "title": "CCS: Read quality",
                    "cpswitch_c_active": False,
                    "data_labels": [
                        {"name": "Reads", "ylab": "Number of reads"},
                        {"name": "Bases", "ylab": "Number of bases"},
                    ],
                },
            ),
        )

    def convert_to_v5(self, data):
        """Convert the v4 format to the new CCS v5 json format"""
        # Initialise the v5 format dictionary
        v5 = {}
        v5["id"] = "ccs_processing"
        attributes = []
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


def split_attributes(attributes):
    """
    Split report attributes into per-group counts (ZMWs, or DS/SS reads with --hd-finder),
    keyed by ID without the group prefix, and the remaining (v8 read) metrics
    """
    groups: dict[str, dict] = {}
    other: dict[str, dict] = {}
    for attr in attributes:
        key = attr["id"].removeprefix("ccs_processing.")
        # CCS before v8 reports ZMW filter reasons without the zmw_ prefix
        if key.startswith("filtered_"):
            key = "zmw_" + key
        for group in GROUPS:
            if key.startswith(group + "_"):
                groups.setdefault(group, {})[key.removeprefix(group + "_")] = attr
                break
        else:
            other[key] = attr
    return groups, other


def parse_PacBio_log(file_content):
    """Parse ccs log file"""
    data = {}
    # This is a local dictionary to store which annotations belong to which
    # result dictionary. This will be used to structure the output, but will
    # not be part of the output itself
    annotations = {}
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
            current_annotation = {}
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
    data = {}

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
