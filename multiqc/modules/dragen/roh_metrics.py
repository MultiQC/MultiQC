'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers ROH metrics data and prepares it for the output report.
It relies on the following official sources:
https://support.illumina.com/help/DRAGEN_Germline_OLH_1000000083701/Content/Source/Informatics/Apps/ROHMetrics_appDRAG.htm
https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/RegionsOfHomozygosity_fDG.htm
https://support-docs.illumina.com/SW/DRAGEN_v41/Content/SW/ROHCaller.htm
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


NAMESPACE = "DRAGEN ROH"


'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The SINGLE_HEADER is used to define common settings for all coverage headers:
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
    "format": "{:,.1f}",  # Value format string - default 1 decimal places
    "cond_formatting_rules": None,  # Rules for conditional formatting table cell values.
    "cond_formatting_colours": None,  # Styles for conditional formatting of table cell values
    "shared_key": None,  # See the link for description
    "modify": None,  # Lambda function to modify values
    "hidden": True,  # Set to True to hide the column in the general table on page load.
    "hidden_own": False,  # For non-general plots in the ROH section.
    "exclude": True,  # True to exclude all headers from the general html table.
}
'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The EXTRA_HEADER contains the logical keys of the SINGLE_HEADER. These are just extensions,
but can be used as if they were a part of the SINGLE_HEADER. Extra configurations are not
set by default. Only the SINGLE_HEADER is used as a basis for a header. But you can always
append an extra config(eg hidden_own, exclude) to the SINGLE_HEADER to overcome this issue.

Hints and rules:
- "hidden" is used for own ROH section if "hidden_own" is not provided.

- "exclude" and "exclude_own" may be useful in cases where html size matters.
  Eg. many samples but only few metrics of interest.

- "order_priority" can be used to arrange table's columns. Only int and float are valid.
'''
EXTRA_HEADER = {
    "hidden_own": None,  # For non-general plots in the ROH section.
    "exclude": None,  # True to exclude metric from the general html table.
    "exclude_own": None,  # True to exclude from own ROH section html table.
    "order_priority": None,  # Used to specify columns' order in all tables.
}

'''"""""""""""""""""""""""""""""""""""""""""""""""
Officially defined file format could not be found.
According to examined real data, there are 4 columns in total:
Section    Void-column    Metric    Value

The only section found in data was 'VARIANT CALLER'
It is unknown if it is the only valid section or others already exist/will be added.
The second column does not contain anything.
It is unknown if it is the final constant structure or it was reserved for future use.

Due to this uncertainty the following dict structure was chosen to store parsed ROH data:
{sample: {section: {region_or_sample: {metric: value}}}}

Might look overcomplicated, but it conforms to the standard.

If you wish to add new sections, take a look at the make_consistent_section.
'''

## Sections:
VARIANT_CALLER = "VARIANT CALLER"

## RG/Sample aka second column:
EMPTY = ""
WGS = "WGS"  # Not detected/used. Just an example.

## Special case. Token for any section/region/sample.
ANY = object()

"""
Only INCLUDED_SECTIONS will be shown in html table.
The parser supports only the first value (4th column).
Can be used as a quick but unreliable way to control html table.
"""
INCLUDED_SECTIONS = {
    VARIANT_CALLER: [EMPTY],
}

# INCLUDED_SECTIONS = {ANY:[ANY]}  # Include all sections/regions/samples.

# Used to detect and report what is not properly supported by the module.
# Works as INCLUDED_SECTIONS.
SUPPORTED_SECTIONS = {
    VARIANT_CALLER: [EMPTY],
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
'''

STD_METRICS = {
    "number of large roh ( >= 3000000)": {
        "order_priority": 1,
        "exclude": False,
        "title": "LargeROH",
        "description": "The number of large ROH detected by the small variant caller.",
        "scale": "Oranges",
        "colour": "255, 255, 0",
        "format": "{:,.0f}",
    },
    "percent snvs in large roh ( >= 3000000)": {
        "order_priority": 2,
        "exclude": False,
        "title": "%SNVsInLargeROH",
        "description": "The percentage of SNVs present in large ROH.",
        "colour": "255, 255, 0",
        "format": "{:,.3f}",
    },
}

# The TABLE_CONFIG defines configs for the whole table in own section.
TABLE_CONFIG = {
    "namespace": NAMESPACE,  # Name for grouping. Prepends desc and is in Config Columns modal
    "id": "dragen_roh_metrics_table",  # ID used for the table
    "table_title": "ROH",  # Title of the table. Used in the column config modal
    "save_file": False,  # Whether to save the table data to a file
    "raw_data_fn": None,  # File basename to use for raw data file
    "sortRows": True,  # Whether to sort rows alphabetically
    "only_defined_headers": True,  # Only show columns that are defined in the headers config
    "col1_header": None,  # The header used for the first column with sample names.
    "no_beeswarm": False,  # Force a table to always be plotted (beeswarm by default if many rows)
}


"""   Define the bargraph's settings   """

CATS = OrderedDict()
CATS["number of large roh ( >= 3000000)"] = {
    "name": "Number of large ROH ( >= 3000000)",
    "color": "#00c3a9",
}
CATS["percent snvs in large roh ( >= 3000000)"] = {
    "name": "Percent SNVs in large ROH ( >= 3000000)",
    "color": "#d2004d",
}

# More configs: https://multiqc.info/docs/development/plots/#bar-graphs
BARGRAPH_CONFIG = {
    "id": "dragen_roh_metrics_bargraph",
    "title": "Dragen: ROH Caller",
    "cpswitch": False,  # Just show the 2 metrics.
    "yDecimals": False,  # Looks a little bit better with integers.
    "ylab": "",  # Just to pass the lint test. The graph is pretty self-explanatory.
    "tt_decimals": -1,  # Preserve values: https://api.highcharts.com/class-reference/Highcharts#.NumberFormatterCallbackFunction
    "tt_percentages": False,  # No need for percentages.
    "stacking": None,  # Both metrics are not parts of something larger.
}


class DragenROHMetrics(BaseMultiqcModule):
    """Public members of the DragenROHMetrics module are available to other dragen modules.
    To avoid cluttering up the global dragen space, only one method is defined."""

    def add_roh_metrics(self):
        """The main function of the dragen ROH metrics module.
        Public members of the BaseMultiqcModule and dragen modules defined in
        MultiqcModule are available within it. Returns keys with sample names
        or an empty set if no samples were found or all samples were ignored."""

        roh_data = {}

        # Stores cleaned sample names with references to file handlers.
        # Used later to check for duplicates and inform about found ones.
        all_samples = defaultdict(list)

        # Metric IDs with original consistent metric strings.
        all_metrics = {}

        for file in self.find_log_files("dragen/roh_metrics"):
            out = roh_parser(file)
            if out["success"]:
                self.add_data_source(file, section="stats")
                cleaned_sample = self.clean_s_name(out["sample_name"], file)
                all_samples[cleaned_sample].append(file)
                roh_data[cleaned_sample] = out["data"]  # Add/overwrite the sample.
                all_metrics.update(out["metric_IDs"])

        # Filter to strip out ignored sample names.
        roh_data = self.ignore_samples(roh_data)
        if not roh_data:
            return set()

        check_duplicate_samples(all_samples, log, "dragen/roh_metrics")

        roh_headers = make_headers(all_metrics)

        txt_out_data = flatten_4D_data(
            roh_data,
            metric_IDs=all_metrics,  # Original strings can be used instead of somewhat ugly metric IDs.
            hide_single_section=True,  # True to exclude common single section from metric ID.
        )
        self.write_data_file(txt_out_data, "dragen_roh")  # Write data to file.

        roh_data_table = make_data_for_table(roh_data)  # Prepare data for tables.
        table_data, table_headers = flatten_4D_data(
            roh_data_table,
            headers=roh_headers,  # Adjust headers' metrics to data's structure.
            metric_IDs=all_metrics,  # Original strings can be used instead of somewhat ugly IDs.
            hide_single_section=True,  # True to exclude common single section from metric ID.
        )
        gen_headers = make_gen_stats_headers(table_headers)
        # If gen_headers is empty, then all metrics from data will be included in table.
        if gen_headers:
            self.general_stats_addcols(table_data, gen_headers, namespace=NAMESPACE)

        # Add single own table section.
        own_table_headers = make_own_headers(table_headers)
        self.add_section(
            name="ROH Caller",
            anchor="dragen-roh-caller-table",
            description="Metrics are reported on a per sample level. Regions of homozygosity (ROH) "
            "are detected as part of the small variant caller. Press the `Help` button for details.",
            helptext="The caller detects and outputs the runs of homozygosity from whole genome "
            "calls on autosomal human chromosomes. ROH output allows downstream tools to screen "
            "for and predict consanguinity between the parents of the proband subject.\n\n"
            "* DRAGEN v3.3 - 3.8: Sex chromosomes are ignored.\n\n"
            "* DRAGEN v3.9 - 4.1: Sex chromosomes are ignored unless the sample sex karyotype is XX, "
            "as specified on the command line or determined by the Ploidy Estimator.\n\n"
            "A region is defined as consecutive variant calls on the chromosome with no large gap in "
            "between these variants. In other words, regions are broken by chromosome or by large gaps "
            "with no SNV calls. The gap size is set to 3 Mbases.\n\nThe ROH algorithm runs on the small "
            "variant calls. The algorithm excludes variants with multiallelic sites, indels, complex "
            "variants, non-PASS filtered calls, and homozygous reference sites. The variant calls are "
            "then filtered further using a block list BED, and finally depth filtering is applied "
            "after the block list filter. The default value for the fraction of filtered calls is 0.2, "
            "which filters the calls with the highest 10% and lowest 10% in DP values. "
            "The algorithm then uses the resulting calls to find regions.\n\n"
            "Addition for DRAGEN v3.7 - 4.1: "
            "The ROH algorithm first finds seed regions that contain at least 50 consecutive homozygous "
            "SNV calls with no heterozygous SNV or gaps of 500,000 bases between the variants. The regions "
            "can be extended using a scoring system that functions as follows.\n\n"
            "* Score increases with every additional homozygous variant (0.025) and decreases with a "
            "large penalty (1–0.025) for every heterozygous SNV. This provides some tolerance of "
            "presence of heterozygous SNV in the region.\n\n"
            "* Each region expands on both ends until the regions reach the end of a chromosome, "
            "a gap of 500,000 bases between SNVs occurs, or the score becomes too low (0).\n\n"
            "Overlapping regions are merged into a single region. Regions can be merged across gaps of "
            "500,000 bases between SNVs if a single region would have been called from the beginning of "
            "the first region to the end of the second region without the gap. There is no maximum size "
            "for regions, but regions always end at chromosome boundaries.",
            plot=table.plot(
                table_data, own_table_headers, {config: val for config, val in TABLE_CONFIG.items() if val is not None}
            ),
        )
        bargraph_data = make_bargraph_data(roh_data)
        if bargraph_data:  # Append only if not empty.
            self.add_section(
                name="Large ROH / SNVs in large ROH",
                anchor="dragen-roh-caller-bargraph",
                description="Different perspective for the 'Number of large ROH ( >= 3000000)' "
                "and 'Percent SNVs in large ROH ( >= 3000000)' metrics from the ROH Caller.",
                plot=bargraph.plot(
                    bargraph_data,
                    CATS,
                    BARGRAPH_CONFIG,
                ),
            )
        # Report found info/warnings/errors. You can disable it anytime, if it is not wanted.
        make_log_report(log_data, log, "roh_metrics")

        return roh_data.keys()


def make_data_for_table(DATA):
    """Modify parsed data before storing in tables."""

    # Inner samples can be concatenated with file samples
    # and replaced by empty string in the second column.

    # Sections can be analyzed and restructured.

    # Use INCLUDED_SECTIONS as inclusion checklist.
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
    data_out = defaultdict(dict)
    included_metrics = ["number of large roh ( >= 3000000)", "percent snvs in large roh ( >= 3000000)"]
    for sample in data:
        # If section is "VARIANT CALLER" and second column is ""
        if VARIANT_CALLER in data[sample] and EMPTY in data[sample][VARIANT_CALLER]:
            metrics = data[sample][VARIANT_CALLER][EMPTY]
            for metric in included_metrics:
                # Crash prevention. Though very unlikely that metric is not present.
                if metric in metrics:
                    data_out[sample][metric] = metrics[metric]
    return data_out


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The log_data is the container for found info, warnings and errors.
Collected data is stored in groups to make the output log file more readable.

warnings:
- invalid_file_names, which do not conform to the:
    <output-file-prefix>.roh_metrics.csv

debug/info:
- invalid_file_lines, which do not conform to:
  <section>,<RG/Sample/Empty>,<metric>,<value>
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
    """Isolation for the parsing codeblock. Returns parse_roh_metrics_file closure."""

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
    # <output-file-prefix>.roh_metrics.csv
    FILE_RGX = re.compile("(.+)\.roh_metrics\.csv$")

    # No official general line structure could be found.
    # The following regex is used to check input lines.
    LINE_RGX = re.compile("^([^,]+),([^,]*),([^,]+),([^,]+)$")

    def parse_roh_metrics_file(file_handler):
        """
        Parser for *.roh_metrics.csv files.

        Input:  file_handler with necessary info - file name/content/root.

        Output: {"success": False} if file name test failed. Otherwise:
                {"success": False/True,
                 "sample_name": <output-file-prefix>,
                 "data": {section: {RG/Sample: {metric: value}}},
                 "metric_IDs": {metric_ID: original_metric}
                }, with success False if all file's lines are invalid.

        File examples:

        sample_1.roh_metrics.csv
        VARIANT CALLER,,Percent SNVs in large ROH ( >= 3000000),0.000
        VARIANT CALLER,,Number of large ROH ( >= 3000000),0

        sample_2.roh_metrics.csv
        VARIANT CALLER,,Percent SNVs in large ROH ( >= 3000000),0.082
        VARIANT CALLER,,Number of large ROH ( >= 3000000),1
        """
        file, root = file_handler["fn"], file_handler["root"]

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
                section, region_sample, metric, value = line_match.groups()

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

            # Check the extracted value. Shall be float or int.
            try:
                value = float(value)
            except ValueError:
                try:
                    value = int(value)
                except ValueError:
                    log_data["unusual_values"][root][file][metric] = value

            metric_IDs[metric_id] = consistent_metric
            data[consistent_section][region_sample][metric_id] = value

        # Return results of parsing the input file.
        return {
            "success": success,
            "sample_name": sample,
            "data": data,
            "metric_IDs": metric_IDs,
        }

    # Return the parser closure.
    return parse_roh_metrics_file


roh_parser = create_parser()
