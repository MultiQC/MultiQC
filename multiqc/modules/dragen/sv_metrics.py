'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers SV metrics data and prepares it for the output report.
It relies on the following official sources (unfortunately uninformative ones, nothing else could be found):
https://support-docs.illumina.com/SW/DRAGEN_v41/Content/SW/DRAGEN/MantaStatisticsOutput_fDG.htm
https://support.illumina.com/help/DRAGEN_Germline_OLH_1000000083701/Content/Source/Informatics/Apps/OutputFiles_swBS_appDRAGGP.htm
'''

import logging
import re
from collections import OrderedDict, defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table
from multiqc.utils.util_functions import write_data_file

from .utils import check_duplicate_samples, flatten_4D_data, make_gen_stats_headers, make_log_report, make_own_headers

# Initialise the logger.
log = logging.getLogger(__name__)


NAMESPACE = "DRAGEN SV"


'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The SINGLE_HEADER is used to define the common settings for all coverage headers:
https://multiqc.info/docs/development/plots/#creating-a-table
'''
SINGLE_HEADER = {
    "namespace": NAMESPACE,  # Name for grouping. Prepends desc and is in Config Columns modal
    "title": None,  # Short title, table column title
    "description": None,  # Longer description, goes in mouse hover text
    "max": None,  # Minimum value in range, for bar / colour coding
    "min": 0,  # Maximum value in range, for bar / colour coding
    "ceiling": None,  # Maximum value for automatic bar limit
    "floor": None,  # Minimum value for automatic bar limit
    "minRange": None,  # Minimum range for automatic bar
    "scale": "GnBu",  # Colour scale for colour coding. False to disable.
    "bgcols": None,  # Dict with values: background colours for categorical data.
    "colour": "0, 0, 255",  # Colour for column grouping
    "suffix": None,  # Suffix for value (eg. "%")
    "format": "{:,.2f}",  # Value format string - Multiqc's default 1 decimal places
    "cond_formatting_rules": None,  # Rules for conditional formatting table cell values.
    "cond_formatting_colours": None,  # Styles for conditional formatting of table cell values
    "shared_key": None,  # See the link for description
    "modify": None,  # Lambda function to modify values
    "hidden": True,  # Set to True to hide the column in the general table on page load.
    "hidden_own": True,  # For non-general plots in the SV section.
    "exclude": True,  # True to exclude all headers from the general html table.
}
'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The EXTRA_HEADER contains the logical keys of the SINGLE_HEADER. These are just extensions,
but can be used as if they were a part of the SINGLE_HEADER. Extra configurations are not
set by default. Only the SINGLE_HEADER is used as a basis for a header. But you can always
append an extra config(eg hidden_own, exclude) to the SINGLE_HEADER to overcome this issue.

Hints and rules:
- "hidden" is used for own SV section if "hidden_own" is not provided.

- "exclude" and "exclude_own" may be useful in cases where html size matters.
  Eg. many samples but only few metrics of interest.

- "order_priority" can be used to arrange table's columns. Only int and float are valid.
'''
EXTRA_HEADER = {
    "hidden_own": None,  # For non-general plots in the SV section.
    "exclude": None,  # True to exclude metric from the general html table.
    "exclude_own": None,  # True to exclude from own SV section html table.
    "order_priority": None,  # Used to specify columns' order in all tables.
}

'''"""""""""""""""""""""""""""""""""""""""""""""""
Officially defined file format could not be found.
According to examined real data, there are 5 columns in total:
Section    RG/Sample    Metric    Value1    Value2

The following structure was chosen to store parsed SV data:
{sample: {section: {RG/Sample: {metric: value}}}}

Might look overcomplicated, but it conforms to the standard.

If you wish to add new sections, take a look at the make_consistent_section.
'''

## Sections:
SV_SUMMARY = "SV SUMMARY"

## RG/Sample aka Second column:
EMPTY = ""
WGS = "WGS"  # Not detected/used. Just an example.

## Special case. Token for any section/region/sample.
ANY = object()

"""
Only INCLUDED_SECTIONS will be shown in html table.
Can be used as a quick but unreliable way to control html table.
"""
INCLUDED_SECTIONS = {
    SV_SUMMARY: [ANY],
}

# INCLUDED_SECTIONS = {ANY:[ANY]}  # Include all sections/regions/samples.

# Used to detect and report what is not properly supported by the module.
# Works as INCLUDED_SECTIONS.
SUPPORTED_SECTIONS = {
    SV_SUMMARY: [ANY],
}

'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The STD_METRICS is the container for abstract representations of the standard metrics.
They are designed to be seen as simple identifiers of the well known metrics.
You can set any real/virtual configuration defined in the SINGLE_HEADER.

Some rules:

- "modify" must conform to the following general structure:
  lambda x: (arbitrary expression) if isinstance(x, str) else (math expression)
  The x is either number or string. The latter is some value, which could not be
  converted to int or float. NA for instance. Such values can be modified in the
  (arbitrary expression) if needed. If not then just return the input.

- If you wish to add new metrics, take a look at the make_consistent_metric.

Note:
  file's sample names are concatenated with internal sample names from the 2nd column.
  That makes them incompatible with sample names produced by other dragen modules,
  so icluding headers into gen stats table would result in a separate row for each concatenation.

  There are at least 2 other ways to handle metrics:
  1) Ignore internal sample if it is a part of the file's sample name. That may result
     in the loss of information because internal sample can have additional info.
     Besides it solves the problem only for the case sample.
  2) Concatenate with metrics. This would make the files' samples consistent across modules,
     but make exclusive additional column for each internal sample.
  Not sure if there is a way to include internal samples into the gen stats table without
  making it ugly and make sample names compatible with other modules.
'''

V2 = " pct"  # Used to create a unique ID if metric has two values.

STD_METRICS = {
    "total number of structural variants (pass)": {
        "order_priority": 1,
        "hidden_own": False,
        "title": "SVs total",
        "description": "Total number of structural variants (PASS).",
        "format": "{:,.0f}",
    },
    "number of deletions (pass)": {
        "order_priority": 2,
        "hidden_own": False,
        "title": "Deletions",
        "description": "Number of deletions (PASS).",
        "scale": "Oranges",
        "format": "{:,.0f}",
    },
    "number of deletions (pass)"
    + V2: {
        "order_priority": 3,
        "title": "Deletions%",
        "description": "Percentage of deletions (PASS).",
        "scale": "Blues",
    },
    "number of insertions (pass)": {
        "order_priority": 4,
        "hidden_own": False,
        "title": "Insertions",
        "description": "Number of insertions (PASS).",
        "scale": "Greens",
        "format": "{:,.0f}",
    },
    "number of insertions (pass)"
    + V2: {
        "order_priority": 5,
        "title": "Insertions%",
        "description": "Percentage of insertions (PASS).",
        "scale": "Oranges",
    },
    "number of duplications (pass)": {
        "order_priority": 6,
        "hidden_own": False,
        "title": "Duplications",
        "description": "Number of duplications (PASS).",
        "scale": "Blues",
        "format": "{:,.0f}",
    },
    "number of duplications (pass)"
    + V2: {
        "order_priority": 7,
        "title": "Duplications%",
        "description": "Percentage of duplications (PASS).",
        "scale": "Greens",
    },
    "number of breakend pairs (pass)": {
        "order_priority": 8,
        "hidden_own": False,
        "title": "BreakendPairs",
        "description": "Number of breakend pairs (PASS).",
        "scale": "Oranges",
        "format": "{:,.0f}",
    },
    "number of breakend pairs (pass)"
    + V2: {
        "order_priority": 9,
        "title": "BreakendPairs%",
        "description": "Percentage of breakend pairs (PASS).",
        "scale": "Blues",
    },
}

# The TABLE_CONFIG defines configs for the whole table in own section.
TABLE_CONFIG = {
    "namespace": NAMESPACE,  # Name for grouping. Prepends desc and is in Config Columns modal
    "id": "dragen_sv_metrics_table",  # ID used for the table
    "table_title": "SV Metrics",  # Title of the table. Used in the column config modal
    "save_file": False,  # Whether to save the table data to a file
    "raw_data_fn": None,  # File basename to use for raw data file
    "sortRows": True,  # Whether to sort rows alphabetically
    "only_defined_headers": True,  # Only show columns that are defined in the headers config
    "col1_header": None,  # The header used for the first column with sample names.
    "no_beeswarm": False,  # Force a table to always be plotted (beeswarm by default if many rows)
}

"""   Define the bargraph's settings   """

CATS = [OrderedDict(), OrderedDict()]
CATS[0]["total number of structural variants (pass)"] = {
    "name": "Total number of structural variants (PASS)",
    "color": "#00c3a9",
}
CATS[1]["number of deletions (pass)"] = {
    "name": "Number of deletions (PASS)",
    "color": "#ffd500",
}
CATS[1]["number of insertions (pass)"] = {
    "name": "Number of insertions (PASS)",
    "color": "#00c4ff",
}
CATS[1]["number of duplications (pass)"] = {
    "name": "Number of duplications (PASS)",
    "color": "#ff004d",
}
CATS[1]["number of breakend pairs (pass)"] = {
    "name": "Number of breakend pairs (PASS)",
    "color": "#00c3a9",
}


# More configs: https://multiqc.info/docs/development/plots/#bar-graphs
BARGRAPH_CONFIG = {
    "id": "sv_metrics_bargraph",
    "title": "Dragen: SV Metrics",
    "tt_decimals": -1,  # Preserve values: https://api.highcharts.com/class-reference/Highcharts#.NumberFormatterCallbackFunction
    "ylab": "",  # Just to pass lint test. 'ylab' is defined in data_labels.
    "data_labels": [
        {
            "name": "Total SVs",
            "ylab": "Counts",
        },
        {
            "name": "SVs Breakdown",
            "ylab": "Counts",
        },
    ],
}


class DragenSVMetrics(BaseMultiqcModule):
    """Public members of the DragenSVMetrics module are available to other dragen modules.
    To avoid cluttering up the global dragen space, only one method is defined."""

    def add_sv_metrics(self):
        """The main function of the dragen SV metrics module.
        Public members of the BaseMultiqcModule and dragen modules defined in
        MultiqcModule are available within it. Returns keys with sample names
        or an empty set if no samples were found or all samples were ignored."""

        sv_data = {}

        # Stores cleaned sample names with references to file handlers.
        # Used later to check for duplicates and inform about found ones.
        all_samples = defaultdict(list)

        # Metric IDs with original consistent metric strings.
        all_metrics = {}

        for file in self.find_log_files("dragen/sv_metrics"):
            out = sv_parser(file)
            if out["success"]:
                self.add_data_source(file, section="stats")
                cleaned_sample = self.clean_s_name(out["sample_name"], file)
                all_samples[cleaned_sample].append(file)
                sv_data[cleaned_sample] = out["data"]  # Add/overwrite the sample.
                all_metrics.update(out["metric_IDs"])

        # Filter to strip out ignored sample names:
        sv_data = self.ignore_samples(sv_data)
        if not sv_data:
            return set()

        check_duplicate_samples(all_samples, log, "dragen/sv_metrics")

        sv_headers = make_headers(all_metrics)

        sv_data_modified = restructure_data(sv_data)

        txt_out_data = flatten_4D_data(
            sv_data_modified,
            metric_IDs=all_metrics,  # Original strings can be used instead of somewhat ugly metric IDs.
            hide_single_section=True,  # True to exclude common single section from metric ID.
        )
        self.write_data_file(txt_out_data, "dragen_sv")  # Write data to file.

        sv_data_table = make_data_for_table(sv_data_modified)  # Prepare data for tables.
        table_data, table_headers = flatten_4D_data(
            sv_data_table,
            headers=sv_headers,  # Adjust headers' metrics to data's structure.
            metric_IDs=all_metrics,  # Original strings can be used instead of somewhat ugly IDs.
            hide_single_section=True,  # True to exclude common single section from metric ID.
        )
        gen_headers = make_gen_stats_headers(table_headers)
        # If gen_headers is empty, then all metrics from data will be included in table.
        if gen_headers:
            self.general_stats_addcols(table_data, gen_headers, namespace=NAMESPACE)

        own_table_headers = make_own_headers(table_headers)
        self.add_section(
            name="SV Metrics",
            anchor="dragen-sv-metrics-table",
            description="Structural Variants (SV) metrics from &#60;output-file-prefix&#62;.sv_metrics.csv files.",
            plot=table.plot(
                table_data, own_table_headers, {config: val for config, val in TABLE_CONFIG.items() if val is not None}
            ),
        )
        bargraph_data = make_bargraph_data(sv_data_modified)
        if bargraph_data:  # Append only if not empty.
            self.add_section(
                name="Structural Variants (SV)",
                anchor="dragen-sv-metrics-bargraph",
                description="Total SVs with Deletions / Insertions / Duplications / Breakend pairs.",
                plot=bargraph.plot(
                    bargraph_data,
                    CATS,
                    BARGRAPH_CONFIG,
                ),
            )
        # Report found info/warnings/errors. You can disable it anytime, if it is not wanted.
        make_log_report(log_data, log, "sv_metrics")

        return sv_data.keys()


def restructure_data(data):
    """Modify data to make it better looking in tables. The output must be 4D."""
    data_out = defaultdict(lambda: defaultdict(dict))
    for sample in data:
        for section in data[sample]:
            # SV_SUMMARY has samples in the second column.
            # These will be concatenated with file's sample name.
            if section == SV_SUMMARY:
                for inner_sample in data[sample][section]:
                    if inner_sample:
                        new_sample = sample + "_" + inner_sample
                    else:
                        new_sample = sample
                    # store data.
                    data_out[new_sample][section][EMPTY] = {
                        metric: value for metric, value in data[sample][section][inner_sample].items()
                    }
            # Other sections (if any exist) are just copied.
            else:
                data_out[sample][section] = {
                    rg_sample: {metric: value for metric, value in metric_values.items()}
                    for rg_sample, metric_values in data[sample][section].items()
                }
    return data_out


def make_data_for_table(DATA):
    """Use INCLUDED_SECTIONS as inclusion checklist."""
    return {
        sample: {
            section: {
                region: {metric: value for metric, value in data.items()}
                for region, data in regions.items()
                if (ANY in INCLUDED_SECTIONS and ANY in INCLUDED_SECTIONS[ANY])
                or (
                    section in INCLUDED_SECTIONS
                    and (region in INCLUDED_SECTIONS[section] or ANY in INCLUDED_SECTIONS[section])
                )
            }
            for section, regions in sections.items()
            if ANY in INCLUDED_SECTIONS or section in INCLUDED_SECTIONS
        }
        for sample, sections in DATA.items()
    }


def make_bargraph_data(data):
    """Extract metrics of interest from restructured data."""
    total_SVs = defaultdict(dict)
    components = defaultdict(dict)
    component_metrics = [
        "number of deletions (pass)",
        "number of insertions (pass)",
        "number of duplications (pass)",
        "number of breakend pairs (pass)",
    ]
    for sample in data:
        for section in data[sample]:
            if section == SV_SUMMARY and EMPTY in data[sample][section]:
                metrics = data[sample][section][EMPTY]
                for metric in metrics:
                    if metric == "total number of structural variants (pass)":
                        total_SVs[sample][metric] = metrics[metric]
                    elif metric in component_metrics:
                        components[sample][metric] = metrics[metric]
    return [total_SVs, components]


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The log_data is the container for found info, warnings and errors.
Collected data is stored in groups to make the output log file more readable.

warnings:
- invalid_file_names, which do not conform to the:
    <output-file-prefix>.sv_metrics.csv

debug/info:
- invalid_file_lines, which do not conform to:
  <section>,<RG/Sample/Empty>,<metric>,<value> or
  <section>,<RG/Sample/Empty>,<metric>,<value1>,<value2>
- unknown_sections are not present in the SUPPORTED_SECTIONS.
- unknown_second_columns are not present in the SUPPORTED_SECTIONS.
- unknown_metrics are not present in the STD_METRICS.
- unusual_values are those except for int and float.
'''
log_data = {
    "invalid_file_names": defaultdict(list),
    "invalid_file_lines": defaultdict(lambda: defaultdict(list)),
    "unknown_sections": set(),
    "unknown_second_columns": defaultdict(set),
    "unknown_metrics": [],
    "unusual_values": defaultdict(lambda: defaultdict(dict)),
}


def make_headers(metric_IDs):
    """Build headers from metric_IDs."""
    headers = {}
    for metric_id in metric_IDs:
        original_metric = metric_IDs[metric_id]
        configs = {config: val for config, val in SINGLE_HEADER.items() if val is not None}
        if metric_id in STD_METRICS:
            configs.update(STD_METRICS[metric_id])
        else:
            log_data["unknown_metrics"].append(original_metric)
            configs["title"] = original_metric
            configs["description"] = original_metric
        headers[metric_id] = configs
    return headers


def create_parser():
    """Isolation for the parsing codeblock. Returns parse_sv_metrics_file closure."""

    def make_consistent_metric(metric):
        """Tries to guarantee consistency to a certain degree.
        Metric-specific peculiarities can be handled here."""
        # metric will be stored in .txt output.
        metric = re.sub("\s+", " ", metric).strip()
        # metric_id will be used for internal storage.
        metric_id = metric.lower()
        return metric, metric_id

    def make_consistent_section(section):
        """Tries to guarantee consistency for section string."""
        return re.sub("\s+", " ", section).strip().upper()

    # The official/accepted structure of the input file:
    # <output-file-prefix>.sv_metrics.csv
    FILE_RGX = re.compile("(.+)\.sv_metrics\.csv$")

    # No official general line structure could be found.
    # The following regex is used to check input lines.
    LINE_RGX = re.compile("^([^,]+),([^,]*),([^,]+),([^,]+)(,[^,]+)?$")

    def parse_sv_metrics_file(file_handler):
        """
        Parser for *.sv_metrics.csv files.

        Input:  file_handler with necessary info - file name/content/root.

        Output: {"success": False} if file name test failed. Otherwise:
                {"success": False/True,
                 "sample_name": <output-file-prefix>,
                 "data": {section: {RG/Sample: {metric: value}}},
                 "metric_IDs": {metric_ID: original_metric}
                }, with success False if all file's lines are invalid.

        File example:

        TTG0_037_dragen.sv_metrics.csv
        SV SUMMARY,TTG0_037,Total number of structural variants (PASS),15972
        SV SUMMARY,TTG0_037,Number of deletions (PASS),6211,38.87
        SV SUMMARY,TTG0_037,Number of insertions (PASS),8935,55.94
        SV SUMMARY,TTG0_037,Number of duplications (PASS),67,0.42
        SV SUMMARY,TTG0_037,Number of breakend pairs (PASS),759,4.75
        SV SUMMARY,TTG0_038,Total number of structural variants (PASS),13289
        SV SUMMARY,TTG0_038,Number of deletions (PASS),5927,44.60
        SV SUMMARY,TTG0_038,Number of insertions (PASS),5691,42.82
        SV SUMMARY,TTG0_038,Number of duplications (PASS),39,0.29
        SV SUMMARY,TTG0_038,Number of breakend pairs (PASS),1632,12.28
        """
        file, root = file_handler["fn"], file_handler["root"]

        # Skip duplicates.
        if file == "diploidSV.sv_metrics.csv":
            return {"success": False}

        file_match = FILE_RGX.search(file)
        if not file_match:
            log_data["invalid_file_names"][root].append(file)
            return {"success": False}

        sample = file_match.group(1)

        success = False
        data = defaultdict(lambda: defaultdict(dict))
        metric_IDs = {}

        for line in file_handler["f"].splitlines():
            # Check the general line structure.
            line_match = LINE_RGX.search(line)

            # If line is fine then extract the necessary fields.
            if line_match:
                success = True
                section, region_sample, metric, value1, value2 = line_match.groups()
                if value2:
                    value2 = value2[1:]  # Get rid of comma.

            # Otherwise check if line is empty. If not then report it. Go to the next line.
            else:
                if not re.search("^\s*$", line):
                    log_data["invalid_file_lines"][root][file].append(line)
                continue

            # Modify extracted parts to improve consistency.
            consistent_section = make_consistent_section(section)
            consistent_metric, metric_id = make_consistent_metric(metric)

            # Check support for sections/regions.
            if not (ANY in SUPPORTED_SECTIONS or consistent_section in SUPPORTED_SECTIONS):
                log_data["unknown_sections"].add(section)
            if not (
                ANY in SUPPORTED_SECTIONS
                and ANY in SUPPORTED_SECTIONS[ANY]
                or (
                    consistent_section in SUPPORTED_SECTIONS
                    and (
                        region_sample in SUPPORTED_SECTIONS[consistent_section]
                        or ANY in SUPPORTED_SECTIONS[consistent_section]
                    )
                )
            ):
                log_data["unknown_second_columns"][section].add(region_sample)

            # Check the extracted values. Must be int or float.
            try:
                value1 = int(value1)
            except ValueError:
                try:
                    value1 = float(value1)
                except ValueError:
                    log_data["unusual_values"][root][file][metric] = value1

            if value2:
                try:
                    value2 = float(value2)
                except ValueError:
                    try:
                        value2 = int(value2)
                    except ValueError:
                        log_data["unusual_values"][root][file][metric + " (second value)"] = value2

            metric_IDs[metric_id] = consistent_metric
            data[consistent_section][region_sample][metric_id] = value1
            if value2 is not None:
                metric_IDs[metric_id + V2] = consistent_metric + V2
                data[consistent_section][region_sample][metric_id + V2] = value2

        # Return results of parsing the input file.
        return {
            "success": success,
            "sample_name": sample,
            "data": data,
            "metric_IDs": metric_IDs,
        }

    # Return the parser closure.
    return parse_sv_metrics_file


sv_parser = create_parser()
